#!/bin/bash

set -e

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Detect operating system and architecture
detect_platform() {
    OS=$(uname -s)
    ARCH=$(uname -m)

    case "$OS" in
        Darwin*)
            PLATFORM="macos"
            if [[ "$ARCH" == "arm64" ]]; then
                echo "Detected: macOS on Apple Silicon (M1/M2)"
                COMPILER_PACKAGES="clang_osx-arm64"
            else
                echo "Detected: macOS on Intel"
                COMPILER_PACKAGES="clang_osx-64"
            fi
            ;;
        Linux*)
            PLATFORM="linux"
            echo "Detected: Linux"
            COMPILER_PACKAGES="gcc=11.3.0"
            ;;
        CYGWIN*|MINGW*|MSYS*)
            PLATFORM="windows"
            echo "Detected: Windows"
            COMPILER_PACKAGES="gcc=11.3.0"
            ;;
        *)
            echo "Unsupported operating system: $OS"
            exit 1
            ;;
    esac
}

# Check if conda is installed
if ! command_exists conda; then
    echo "Conda is not installed. Please install Conda first."
    echo "Visit: https://docs.conda.io/projects/conda/en/latest/user-guide/install/"
    exit 1
fi

detect_platform

# Set up conda env: create if missing, otherwise sync to environment.yml
echo "Setting up conda environment..."
if conda env list | grep -qE "^bcar-env\s"; then
    echo "Environment 'bcar-env' exists; syncing to environment.yml..."
    conda env update -n bcar-env -f environment.yml --prune
else
    echo "Creating environment 'bcar-env'..."
    conda env create -f environment.yml
fi

# Make sure the compiler is installed
conda install -n bcar-env -c conda-forge $COMPILER_PACKAGES -y

echo "Activating environment..."
eval "$(conda shell.bash hook)"
conda activate bcar-env

# Set platform-specific compilers
case "$PLATFORM" in
    macos)
        export CC=clang
        echo "Using Clang compiler for macOS"
        ;;
    linux|windows)
        export CC=gcc
        echo "Using GCC compiler for Linux/Windows"
        ;;
esac

# Display compiler information
echo "Compiler information:"
echo "CC: $CC ($(which $CC))"
$CC --version | head -1 2>/dev/null || echo "Warning: Could not get CC version"

# Clean and compile
echo "Compiling..."
make clean
if make; then
    echo "Installation complete!"
    echo ""
    echo "To use this environment in the future, run:"
    echo "conda activate bcar-env"
else
    echo "Compilation failed. Please check the error messages above."
    exit 1
fi
