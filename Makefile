##################################################
##         MAIN AZTEKAS PROJECT MAKEFILE         #
##################################################

##################################################
# Global Include
##################################################
# Include Makefiles from subdirectories
include $(AZTEKAS_PATH)/include/Makefile
include $(AZTEKAS_PATH)/src/Makefile

.PHONY: init clean pre-commit commit docs help

##################################################
# Environment Management
##################################################
# Install dependencies and activate the virtual environment with Poetry

## Install dependencies and activate the virtual environment
init: aztekas-simulation
	@echo "Installing dependencies and activating virtual environment..."
	@poetry config virtualenvs.path '${AZTEKAS_PATH}/.venv-dev' --local
	@poetry install --no-root
	@poetry shell

aztekas-simulation:
	@echo ""
	@echo "\033[1;39m#################################\033[0m"
	@echo "\033[1;39m######### AZTEKAS-CODE ##########\033[0m"
	@echo "\033[1;39m#################################\033[0m"
	@echo ""

##################################################
# Cleaning Targets
##################################################
# Clean temporary files and executables

## Clean temporary files and executables
clean: aztekas-simulation
	@echo "Cleaning temporary files..."
	rm -rf obj/*.o $(EXECUTABLE) last_$(EXECUTABLE)

##################################################
# Code Quality and Git
##################################################
# Manage pre-commit hooks and Git operations

## Run pre-commit hooks on all files
pre-commit: aztekas-simulation
	@echo "Running pre-commit hooks..."
	@git add $(AZTEKAS_PATH)
	@pre-commit run

## Perform a commit with Commitizen
commit: pre-commit
	@echo "Performing a commit with Commitizen..."
	@cz commit

##################################################
# Documentation
##################################################
# Generate project documentation

## Generate documentation with Doxygen
docs:
	@echo "Generating documentation with Doxygen..."
	doxygen docs/Doxyfile

##################################################
# Help System
##################################################
# Display available commands and their descriptions

## Display available commands and descriptions
help: aztekas-simulation
	@echo "Available commands:"
	@echo "init:         Install dependencies and activate the virtual environment with Poetry"
	@echo "clean:        Clean temporary files and executables"
	@echo "pre-commit:   Run pre-commit hooks on all files"
	@echo "commit:       Perform a Git commit with Commitizen"
	@echo "docs:         Generate project documentation with Doxygen"
	@echo "help:         Display this help message"
