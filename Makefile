# Makefile en la raíz: mi_proyecto/Makefile
# Solo tiene instrucciones de desarrollo e invoca los Makefiles de subdirectorios

include $(AZTEKAS_PATH)/include/Makefile
include $(AZTEKAS_PATH)/src/Makefile

.PHONY: pre-commit commit docs help

## Realiza un commit rápido (como ejemplo). 
## Puedes personalizar el mensaje o quitar este target si no lo necesitas.
clean:
	rm -rf obj/*.o $(EXECUTABLE)

commit:
	cz commit

## Genera documentación con Doxygen
docs:
	doxygen docs/Doxyfile

## Ejecuta pre-commit (formatos, linters, etc.)
pre-commit:
	git add $(AZTEKAS_PATH)
	pre-commit run 

## Ayuda
help:
	@echo "Comandos disponibles:"
	@echo "  make all (o make build)   -> Compila todo el proyecto"
	@echo "  make clean                -> Limpia archivos intermedios"
	@echo "  make pre-commit           -> Ejecuta hooks pre-commit en todos los archivos"
	@echo "  make commit               -> Realiza git add y git commit con un mensaje genérico"
	@echo "  make docs                 -> Genera la documentación con Doxygen"
	@echo "  make help                 -> Muestra esta ayuda"
