#!/bin/sh

check_venv() {
    # python -c "import sys; print(sys.prefix != sys.base_prefix)"
    if [ "$VIRTUAL_ENV" = "" ]; then
        echo "Must be in a virtual environment to update."
        exit 1
    fi
}

remove_directory() {
    find . -name "$1" -type d \
        -exec echo "removing {}" \; \
        -exec rm -dr {} +
}

remove_file() {
    find . -name "$1" -type f \
        -exec echo "removing {}" \; \
        -exec rm {} +
}

case "$1" in
-c)
    remove_directory "__pycache__"
    remove_directory ".ipynb_checkpoints"
    remove_directory ".mypy_cache"
    remove_directory ".ruff_cache"
    ;;
-f)
    echo "+-------------+"
    echo "| ruff format |"
    echo "+-------------+"
    python -m ruff format .
    ;;
-t)
    echo "+------+"
    echo "| mypy |"
    echo "+------+"
    python3 -m mypy ./src

    echo "+------------+"
    echo "| ruff check |"
    echo "+------------+"
    python3 -m ruff check ./src
    ;;
-u)
    check_venv
    python3 -m pip install --upgrade pip
    python3 -m pip install --upgrade -r requirements.txt
    ;;
-j)
    check_venv
    jupyter lab
    ;;
*)
    echo "The choice are:"
    echo "  > [-c] for cleaning the temporary python file;"
    echo "  > [-f] for formatting the code;"
    echo "  > [-j] lauch Jupyter;"
    echo "  > [-t] for testing the code;"
    echo "  > [-u] for updating the python package."
    exit 1
    ;;
esac
