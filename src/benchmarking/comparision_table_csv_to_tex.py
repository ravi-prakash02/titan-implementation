import pandas as pd

# Read CSV
df = pd.read_csv("comparison_table.csv")

# Reorder + rename columns
df = df[[
    "m",
    "Scheme",
    "Commit_ms",
    "Prove_ms",
    "Verify_ms",
    "Proof_bytes",
]]

df = df.rename(columns={
    "Commit_ms": "Commit (ms)",
    "Prove_ms": "Eval (ms)",
    "Verify_ms": "Verify (ms)",
    "Proof_bytes": "Proof (B)",
})

# Sort by m, then Scheme
df = df.sort_values(by=["m", "Scheme"])

# Begin LaTeX table
lines = []
lines.append(r"\begin{table}[t]")
lines.append(r"\centering")
lines.append(r"\caption{Comparison of polynomial commitment schemes grouped by number of variables $m$.}")
lines.append(r"\label{tab:pcs-comparison}")
lines.append(r"\small")
lines.append(r"\begin{tabular}{clrrrr}")
lines.append(r"\toprule")
lines.append(r"$m$ & Scheme & Commit (ms) & Eval (ms) & Verify (ms) & Proof (B) \\")
lines.append(r"\midrule")

# Generate multirow blocks
for m, group in df.groupby("m"):
    rows = group.to_dict("records")
    k = len(rows)

    for i, row in enumerate(rows):
        if i == 0:
            lines.append(
                rf"\multirow{{{k}}}{{*}}{{{m}}} & "
                rf"{row['Scheme']} & "
                rf"{row['Commit (ms)']:.2f} & "
                rf"{row['Eval (ms)']:.2f} & "
                rf"{row['Verify (ms)']:.2f} & "
                rf"{int(row['Proof (B)'])} \\"
            )
        else:
            lines.append(
                rf" & {row['Scheme']} & "
                rf"{row['Commit (ms)']:.2f} & "
                rf"{row['Eval (ms)']:.2f} & "
                rf"{row['Verify (ms)']:.2f} & "
                rf"{int(row['Proof (B)'])} \\"
            )

    lines.append(r"\addlinespace")

lines.append(r"\bottomrule")
lines.append(r"\end{tabular}")
lines.append(r"\end{table}")

# Write to file
with open("comparison_table.tex", "w") as f:
    f.write("\n".join(lines))

print("LaTeX table with multirow generated successfully.")
