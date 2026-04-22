import matplotlib.pyplot as plt
import pandas as pd

def load_scheme(path):
    df = pd.read_csv(path, header=None)
    return {
        "m": df.iloc[:, 0],
        "commit": df.iloc[:, 1],
        "eval": df.iloc[:, 2],
        "verify": df.iloc[:, 3],
        "proof_size": df.iloc[:, 4],
    }

schemes = {
    "Titan": load_scheme("titan_clean.csv"),
    "Dory": load_scheme("dory_BLS_clean.csv"),
    "WHIR": load_scheme("whir_clean.csv"),
    "Kopis": load_scheme("kopis_clean.csv"),
    "Brakedown": load_scheme("brakedown_clean.csv"),
    "hyrax": load_scheme("hyrax_clean.csv"),
}

def plot_metric(metric, ylabel, filename):
    plt.figure()

    for name, data in schemes.items():
        plt.plot(
            data["m"],
            data[metric],
            marker="o",
            label=name,
        )

    plt.xlabel("Number of variables (m)")
    plt.ylabel(ylabel)
    plt.yscale("log")
    plt.legend()
    plt.grid(True, which="both")
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

plot_metric("commit", "Commit time (ms)", "commit_time.png")
plot_metric("eval", "Eval time (ms)", "eval_time.png")
plot_metric("verify", "Verify time (ms)", "verify_time.png")
plot_metric("proof_size", "Proof size (bytes)", "proof_size.png")


rows = []

for scheme_name, data in schemes.items():
    for i in range(len(data["m"])):
        rows.append({
            "Scheme": scheme_name,
            "m": data["m"].iloc[i],
            "Commit_ms": data["commit"].iloc[i],
            "Eval_ms": data["eval"].iloc[i],
            "Verify_ms": data["verify"].iloc[i],
            "Proof_bytes": data["proof_size"].iloc[i],
        })

full_table = pd.DataFrame(rows)

# Save CSV
full_table.to_csv("comparison_table.csv", index=False)
df = pd.read_csv("comparison_table.csv")

print("All plots and tables generated successfully.")
