from Bio import Entrez, SeqIO
from http.client import IncompleteRead
import pandas as pd
import matplotlib.pyplot as plt
import time
import random

printer = print

class GenBankRetriever:
    def __init__(self, email, api_key):
        Entrez.email = email
        Entrez.api_key = api_key
        self.email = email
        self.api_key = api_key
        self.webenv = None
        self.query_key = None

    def search_by_taxid(self, taxid):
   
        try:
            tax_handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            tax_record = Entrez.read(tax_handle)
            organism = tax_record[0]['ScientificName']
            printer(f"Organism found: {organism}")

            query = f"txid{taxid}[Organism]"
            handle = Entrez.esearch(db="nucleotide", term=query, usehistory="y")
            search_result = Entrez.read(handle)

            self.webenv = search_result["WebEnv"]
            self.query_key = search_result["QueryKey"]
            total_count = int(search_result["Count"])
            return total_count

        except Exception as e:
            printer(f"Error during taxid search: {e}")
            return 0

    def fetch_records(self, max_records, min_len, max_len):
        start = 0
        batch_size = 500
        data = []

        while start < max_records:
            try:
                handle = Entrez.efetch(
                    db="nucleotide",
                    rettype="gb",
                    retmode="text",
                    retstart=start,
                    retmax=min(batch_size, max_records - start),
                    webenv=self.webenv,
                    query_key=self.query_key
                )

                for record in SeqIO.parse(handle, "gb"):
                    length = len(record.seq)
                    if min_len <= length <= max_len:
                        data.append((record.id, length, record.description))

                start += batch_size
                printer(f"Downloaded: {min(start, max_records)}/{max_records}")

                time.sleep(0.5 + random.random() * 0.3)

            except IncompleteRead:
                printer("IncompleteRead encountered, retrying after pause...")
                time.sleep(5)
                continue

        return data


def save_to_csv(data, filename):
    df = pd.DataFrame(data, columns=["Accession", "Length", "Description"])
    df.to_csv(filename, index=False)
    printer(f"Saved CSV to: {filename}")


def create_plot(csv_file, output_image):
    df = pd.read_csv(csv_file).sort_values("Length", ascending=False).head(100)
    plt.figure(figsize=(12, 6))
    plt.plot(df["Accession"], df["Length"], marker='o')
    plt.xticks(rotation=90, fontsize=6)
    plt.xlabel("Accession")
    plt.ylabel("Sequence Length")
    plt.title("Top 100 Longest Sequences")
    plt.tight_layout()
    plt.savefig(output_image)
    printer(f"Plot saved to: {output_image}")


def main():
    printer("Enter your email:")
    email = input()
    printer("Enter your API key:")
    api_key = input()
    printer("Enter TaxID:")
    taxid = input()
    printer("Minimum sequence length:")
    min_length = int(input())
    printer("Maximum sequence length:")
    max_length = int(input())
    printer("Maximum number of records to fetch:")
    max_records = int(input())

    retriever = GenBankRetriever(email, api_key)
    total = retriever.search_by_taxid(taxid)
    printer(f"Total records found: {total}")

    data = retriever.fetch_records(max_records, min_length, max_length)
    printer(f"Filtered matches found: {len(data)}")

    csv_filename = f"{taxid}_{min_length}_{max_length}.csv"
    png_filename = f"{taxid}_{min_length}_{max_length}.png"
    save_to_csv(data, csv_filename)
    create_plot(csv_filename, png_filename)

if __name__ == "__main__":
    main()
