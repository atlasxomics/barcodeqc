import os
import pandas as pd
from pathlib import Path

for f in os.listdir("merLists"):
    print(f)
    f = Path(f"merLists/{f}")
    name = f.name

    data_types = {
        'Plate': str,
        'plateRow': str,
        'plateCol': int,
        'DNA_name': str,
        '8mer': str,
        'spatialRow': 'Int32',
        'spatialCol': 'Int32'
    }
    df = pd.read_csv(f, sep="\t", header=0, dtype=data_types)
    df = df.dropna()
    df = df.drop(["plateRow", "plateCol", "DNA_name"], axis=1)

    df = df.rename(
        columns={"8mer": "sequence", "spatialCol": "col", "spatialRow": "row"}
    )

    df_a = df[df['Plate'].str.contains("BCA")]
    df_b = df[df['Plate'].str.contains("BCB")]

    df_a = df_a.drop(["Plate"], axis=1)
    df_b = df_b.drop(["Plate"], axis=1)

    df_a["col"] = df_a["col"].astype(int)
    df_a["row"] = df_a["row"].astype(int)

    df_b["col"] = df_b["col"].astype(int)
    df_b["row"] = df_b["row"].astype(int)

    df_a.to_csv(f"barcode_files/barcode_A/{name}", index=False)
    df_b.to_csv(f"barcode_files/barcode_B/{name}", index=False)
