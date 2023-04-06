from db import get_mongo_client
import numpy as np
import click
@click.command()
@click.argument("pid")
def main(pid):
    pid = np.nan
    db_client = get_mongo_client()
    db = db_client['ccgp_dev']
    metadata = db["sample_metadata"]

    docs = list(metadata.find({"ccgp-project-id": pid}))
    for d in docs:
        samp = d["*sample_name"]
        l = list(metadata.find({"*sample_name": samp}))
        print(f"Found {len(l)} samples in: {samp}")
    print(f"Would delete {len(docs)} samples from project: {pid}")
    if click.confirm('Do you want to continue?', default=True):
        d = metadata.delete_many({"ccgp-project-id": pid}) 
        print(f"Deleted the following {len(list(d))} documents:")
        for doc in d:
            print(doc)
    else:
        return


if __name__ == "__main__":
    main()