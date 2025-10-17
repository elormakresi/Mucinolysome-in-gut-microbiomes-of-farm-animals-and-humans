import json
import pandas as pd
# Step 1: Read file paths from files.txt
with open("./AF3/selected_json/json_file_paths.txt", "r") as f:
    file_paths = [line.strip() for line in f if line.strip()]
# Step 2: Load each JSON and store as row
records = []
for path in file_paths:
    try:
        with open(path, "r") as jf:
            data = json.load(jf)
            if isinstance(data, dict):
                data["source_file"] = path  # Optionally track source
                records.append(data)
            else:
                print(f"Skipped (not dict): {path}")
    except Exception as e:
        print(f"Failed to read {path}: {e}")
# Step 3: Convert to table
df = pd.DataFrame(records)
# Step 4: Show or save
print(df.head())  # Show sample
#df.to_csv("combined_output.csv", index=False)  # Optional: Save as CSV
df.to_csv("other_ppi_scores_mutput.tsv", sep="\t", index=False)  # Save as TSV
