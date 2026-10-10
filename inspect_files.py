import requests
import zipfile
import io

urls = {
    "profiles": "https://zenodo.org/record/14261087/files/SMG_profiles.zip?download=1",
    "metadata": "https://zenodo.org/record/14261087/files/Metadata_All.zip?download=1"
}

def inspect(url, file_path):
    print(f"\n--- Inspecting {file_path} from {url} ---")
    try:
        r = requests.get(url)
        with zipfile.ZipFile(io.BytesIO(r.content)) as z:
            with z.open(file_path) as f:
                for i in range(5):
                    print(f.readline().decode('utf-8', errors='replace').strip())
    except Exception as e:
        print(f"Error: {e}")

inspect(urls["profiles"], "SMG/Taxonomy/mOTUs.txt")
inspect(urls["metadata"], "Metadata/metadata.smg.txt")
