# Makefile en la raíz: mi_proyecto/Makefile
# Solo tiene instrucciones de desarrollo e invoca los Makefiles de subdirectorios

.PHONY: all build clean pre-commit commit docs help

##################################################
# PATH TO AZTEKAS
##################################################

AZTPATH = $(AZTEKAS_PATH)


## Este target compila todo el proyecto llamando a los Makefiles de subdirectorios
all: build

## Construcción completa del proyecto
build:
	@$(MAKE) -C $(AZTPATH)/include
	@$(MAKE) -C $(AZTPATH)/include/aztekas
	@$(MAKE) -C $(AZTPATH)/src

## Limpieza de todos los binarios/object files intermedios
clean:
	@$(MAKE) -C $(AZTPATH)/include clean
	@$(MAKE) -C $(AZTPATH)/include/aztekas clean
	@$(MAKE) -C $(AZTPATH)/src clean

## Ejecuta pre-commit (formatos, linters, etc.)
pre-commit:
	git add $(AZTPATH)
	pre-commit run 

## Realiza un commit rápido (como ejemplo). 
## Puedes personalizar el mensaje o quitar este target si no lo necesitas.
commit:
	cz commit

## Genera documentación con Doxygen
docs:
	doxygen docs/Doxyfile

## Ayuda
help:
	@echo "Comandos disponibles:"
	@echo "  make all (o make build)   -> Compila todo el proyecto"
	@echo "  make clean                -> Limpia archivos intermedios"
	@echo "  make pre-commit           -> Ejecuta hooks pre-commit en todos los archivos"
	@echo "  make commit               -> Realiza git add y git commit con un mensaje genérico"
	@echo "  make docs                 -> Genera la documentación con Doxygen"
	@echo "  make help                 -> Muestra esta ayuda"
