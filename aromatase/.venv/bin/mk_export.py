#!/opt/msc/msc_bioinformatics/aromatase/.venv/bin/python3
import sys
from meeko.cli.mk_export import main
if __name__ == '__main__':
    if sys.argv[0].endswith('.exe'):
        sys.argv[0] = sys.argv[0][:-4]
    sys.exit(main())
