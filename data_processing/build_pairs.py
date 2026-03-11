import argparse
from pathlib import Path
import numpy as np
import pandas as pd

from soundvectors import SoundVectors
from pyclts import CLTS


INVALID_SEGMENTS = {"ː", "0"}

PAIR_SPECS = [
    ("C1_curr", "C2_curr", "12"),
    ("C2_curr", "C3_curr", "23"),
    ("C1_curr", "C3_curr", "13"),
]


def cosine_similarity(v1: np.ndarray, v2: np.ndarray) -> float:
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 == 0 or n2 == 0:
        return np.nan
    return float(np.dot(v1, v2) / (n1 * n2))


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


def normalize_segment(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def normalize_pair_string(s1: str, s2: str) -> str:
    return f"{s1} {s2}"


def normalize_root_string(c1: str, c2: str, c3: str) -> str:
    return f"{c1} {c2} {c3}"


def get_vec_cache(segments, sv):
    vec_cache = {}
    invalid = set()

    for seg in segments:
        seg = normalize_segment(seg)

        if not seg:
            vec_cache[seg] = None
            continue

        if seg in INVALID_SEGMENTS:
            vec_cache[seg] = None
            invalid.add(seg)
            continue

        try:
            vec = sv.get_vec(seg)
            if vec is None:
                vec_cache[seg] = None
                invalid.add(seg)
            else:
                vec_cache[seg] = np.asarray(vec, dtype=float)
        except Exception:
            vec_cache[seg] = None
            invalid.add(seg)

    return vec_cache, invalid


def parse_consonants_column(df: pd.DataFrame) -> pd.DataFrame:
    parts = df["consonants"].str.split(r"\s+", n=1, expand=True)
    if parts.shape[1] != 2:
        raise ValueError("Could not split 'consonants' into exactly two segments.")
    out = df.copy()
    out["C1"] = parts[0]
    out["C2"] = parts[1]
    return out


def build_observed_pairs(df_roots: pd.DataFrame, vec_cache: dict) -> pd.DataFrame:
    rows = []

    for _, row in df_roots.iterrows():
        c1 = row["C1_curr"]
        c2 = row["C2_curr"]
        c3 = row["C3_curr"]

        root = normalize_root_string(c1, c2, c3)

        for left_col, right_col, position in PAIR_SPECS:
            s1 = row[left_col]
            s2 = row[right_col]

            v1 = vec_cache.get(s1)
            v2 = vec_cache.get(s2)

            # Keep only rows where both segments are valid
            if v1 is None or v2 is None:
                continue

            sim = cosine_similarity(v1, v2)

            rows.append(
                {
                    "root": root,
                    "consonants": normalize_pair_string(s1, s2),
                    "position": position,
                    "similarity": sim,
                }
            )

    out = pd.DataFrame(rows, columns=["root", "consonants", "position", "similarity"])
    return out


def build_inventory(df_roots: pd.DataFrame) -> list[str]:
    inventory = pd.unique(df_roots[["C1_curr", "C2_curr", "C3_curr"]].values.ravel("K"))
    inventory = [normalize_segment(x) for x in inventory]
    inventory = [x for x in inventory if x and x not in INVALID_SEGMENTS]
    return list(inventory)


def build_all_possible_pairs_by_position(
    df_roots: pd.DataFrame,
    observed_df: pd.DataFrame,
    vec_cache: dict
) -> pd.DataFrame:
    inventory = build_inventory(df_roots)
    positions = ["12", "23", "13"]

    # Full cartesian product: C1 x C2 x position
    full_index = pd.MultiIndex.from_product(
        [inventory, inventory, positions],
        names=["C1", "C2", "position"]
    )
    all_pairs = full_index.to_frame(index=False)

    # Compute similarity for each pair
    sims = []
    for _, row in all_pairs.iterrows():
        c1 = row["C1"]
        c2 = row["C2"]

        v1 = vec_cache.get(c1)
        v2 = vec_cache.get(c2)

        if v1 is None or v2 is None:
            sims.append(np.nan)
        else:
            sims.append(cosine_similarity(v1, v2))

    all_pairs["similarity"] = sims

    # Count observed frequencies by pair and position
    observed_parsed = parse_consonants_column(observed_df)

    observed_counts = (
        observed_parsed
        .groupby(["C1", "C2", "position"], dropna=False)
        .size()
        .reset_index(name="n")
    )

    out = all_pairs.merge(observed_counts, on=["C1", "C2", "position"], how="left")
    out["n"] = out["n"].fillna(0).astype(int)
    out["occurs"] = (out["n"] > 0).astype(int)
    out["freq"] = out["n"]

    return out[["C1", "C2", "position", "similarity", "n", "occurs", "freq"]]


def main():
    parser = argparse.ArgumentParser(
        description=(
            "From a roots file with C1_curr/C2_curr/C3_curr, build "
            "(1) observed consonant-pair similarities and "
            "(2) all possible consonant-pair combinations by position."
        )
    )
    parser.add_argument("--input", required=True, help="Input CSV/TSV file with C1_curr, C2_curr, C3_curr")
    parser.add_argument("--similarities-out", required=True, help="Output file for observed similarities")
    parser.add_argument("--all-pairs-out", required=True, help="Output file for all pair-position combinations")
    args = parser.parse_args()

    input_file = Path(args.input)
    similarities_out = Path(args.similarities_out)
    all_pairs_out = Path(args.all_pairs_out)

    df = read_table(input_file)

    required_cols = {"C1_curr", "C2_curr", "C3_curr"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    for col in ["C1_curr", "C2_curr", "C3_curr"]:
        df[col] = df[col].map(normalize_segment)

    bipa = CLTS().bipa
    sv = SoundVectors(ts=bipa)

    unique_segments = pd.unique(df[["C1_curr", "C2_curr", "C3_curr"]].values.ravel("K"))
    vec_cache, invalid_segments = get_vec_cache(unique_segments, sv)

    # 1. observed pair similarities
    observed_df = build_observed_pairs(df, vec_cache)
    observed_df.to_csv(similarities_out, sep=";", index=False, encoding="utf-8")

    # 2. all possible pair-position combinations
    all_pairs_df = build_all_possible_pairs_by_position(df, observed_df, vec_cache)
    all_pairs_df.to_csv(all_pairs_out, sep=";", index=False, encoding="utf-8")

    print(f"Wrote {len(observed_df)} rows to {similarities_out}")
    print(f"Wrote {len(all_pairs_df)} rows to {all_pairs_out}")

    if invalid_segments:
        print("\nUnrecognized / skipped segments:")
        for seg in sorted(invalid_segments):
            print(f"  {repr(seg)}")


if __name__ == "__main__":
    main()
