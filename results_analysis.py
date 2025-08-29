import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# Load results
df = pd.read_csv('/wistar/auslander/Timothy/microproteins_timka/hyperparameter_results_negative_num_training.csv')

# Sort if needed
df = df.sort_values(by=['TRAINING'])

# Plot
plt.figure(figsize=(8,6))
sns.regplot(
    data=df,
    x="TRAINING", 
    y="AUROC",
    scatter_kws={"s":60, "color":"blue"},  
    line_kws={"color":"red", "linestyle":"--"},  
    order=2  
)

plt.title("AUROC vs Training Size")
plt.xlabel("Training Size")
plt.ylabel("AUROC")
plt.grid(True)

# Save instead of show
outpath = "/wistar/auslander/Timothy/microproteins_timka/auroc_vs_training.png"
plt.savefig(outpath, dpi=300, bbox_inches="tight")
plt.close()

print(f"Plot saved to {outpath}")