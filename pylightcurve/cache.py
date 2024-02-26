import click


@click.command()
@click.option("--cache-dir", type=click.Path(), default="~/.lightcurve")
def download_database_cache(cache_dir: str):
    """Download the cache from the internet."""
    from pylightcurve.__databases__ import PlcData
    import pathlib

    cache_dir = str(pathlib.Path(cache_dir).expanduser())

    click.echo(f"Downloading cache to {cache_dir}")
    PlcData(
        _reset=True,
        cache_dir=cache_dir,
        ignore_check=False,
    )
    click.echo("Cache downloaded")
    click.echo("Set PYLC_CACHE_DIR to the cache directory.")
    click.echo("Set PYLC_IGNORE_CHECK to True to ignore the database check.")


if __name__ == "__main__":
    download_database_cache()
