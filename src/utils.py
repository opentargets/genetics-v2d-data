from pathlib import Path


def get_project_root() -> Path:
    return Path.cwd()


def get_config_path() -> Path:
    return 'configs' / 'config.yaml'
