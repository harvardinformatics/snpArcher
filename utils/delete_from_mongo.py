from db import get_mongo_client
import click
@click.command()
@click.argument("pid")
def main(pid):
    
    db_client = get_mongo_client()
    db = db_client['ccgp_dev']
    metadata = db["sample_metadata"]

    docs = list(metadata.find({"ccgp-project-id": pid}))

    print(f"Would delete {len(docs)} samples from project: {pid}")
    if click.confirm('Do you want to continue?', default=True):
        d = metadata.delete_many({"ccgp-project-id": pid}) 
        print(f"Deleted the following {len(d)} documents:")
        for doc in d:
            print(doc)
    else:
        return


if __name__ == "__main__":
    main()