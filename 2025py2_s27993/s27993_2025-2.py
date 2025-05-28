from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import time
import random
from http.client import IncompleteRead

def initialize_entrez(email, api_key):
    Entrez.email = email
    Entrez.api_key = api_key
    Entrez.tool = "GenBankExtendedTool"

def search_genbank_by_taxid(taxid):
    try:
        taxonomy_response = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        taxonomy_data = Entrez.read(taxonomy_response)
        organism = taxonomy_data[0]["ScientificName"]
        print(f"Organism: {organism}")
        
        query = f"txid{taxid}[Organism]"
        search_response = Entrez.esearch(db="nucleotide", term=query, usehistory="y")
        search_data = Entrez.read(search_response)
        total_records = int(search_data["Count"])
        
        return total_records, search_data["WebEnv"], search_data["QueryKey"]
    
    except Exception as error:
        print(f"An error occurred during search: {error}")
        return 0, None, None

def fetch_filtered_records(total, webenv, query_key, max_fetch, min_len, max_len):
    filtered = []
    start = 0
    batch_size = 500

    while start < max_fetch:
        try:
            print(f"Fetching records {start + 1} to {min(start + batch_size, max_fetch)}")
            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=min(batch_size, max_fetch - start),
                webenv=webenv,
                query_key=query_key
            )
            for record in SeqIO.parse(handle, "gb"):
                length = len(record.seq)
                if min_len <= length <= max_len:
                    filtered.append((record.id, length, record.description))
            start += batch_size
            time.sleep(0.4 + random.uniform(0, 0.3))
        except IncompleteRead:
            print("Incomplete read occurred. Retrying after a short wait...")
            time.sleep(3)
            continue
    return filtered

def export_to_csv(data, filename):
    df = pd.DataFrame(data, columns=["Accession", "Length", "Description"])
    df.to_csv(filename, index=False)
    print(f"CSV saved: {filename}")
    return df

def generate_plot(df, image_filename):
    df_sorted = df.sort_values("Length", ascending=False).head(100)
    plt.figure(figsize=(12, 6))
    plt.plot(df_sorted["Accession"], df_sorted["Length"], marker='o', linestyle='-')
    plt.xticks(rotation=90, fontsize=6)
    plt.title("Top 100 GenBank Sequences by Length")
    plt.xlabel("Accession Number")
    plt.ylabel("Sequence Length")
    plt.tight_layout()
    plt.savefig(image_filename)
    print(f"Plot saved: {image_filename}")

def main():
    print("Enter your NCBI email:")
    email = input().strip()
    print("Enter your NCBI API key:")
    api_key = input().strip()
    initialize_entrez(email, api_key)

    print("Enter organism's taxonomic ID:")
    taxid = input().strip()

    print("Minimum sequence length:")
    min_length = int(input())
    print("Maximum sequence length:")
    max_length = int(input())
    print("Maximum records to fetch:")
    max_records = int(input())

    total, webenv, query_key = search_genbank_by_taxid(taxid)
    if total == 0 or not webenv or not query_key:
        print("No results found or an error occurred.")
        return

    records = fetch_filtered_records(total, webenv, query_key, max_records, min_length, max_length)
    csv_file = f"{taxid}_{min_length}_{max_length}.csv"
    img_file = f"{taxid}_{min_length}_{max_length}.png"

    df = export_to_csv(records, csv_file)
    generate_plot(df, img_file)

if __name__ == "__main__":
    main()
