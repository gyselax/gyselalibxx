import argparse
import json
from pathlib import Path

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file_in', type=str)
    parser.add_argument('config_file_out', type=str)
    parser.add_argument('relevant_folders', type=Path, nargs='+')
    args = parser.parse_args()

    with open(args.config_file_in, encoding='utf-8') as f:
        config = json.load(f)

    relevant_folders = [r.absolute() for r in args.relevant_folders]

    relevant_configs = [c for c in config if any(r in Path(c['file']).parents for r in relevant_folders)]

    with open(args.config_file_out, 'w', encoding='utf-8') as f:
        config = json.dump(relevant_configs, f, indent=2)
