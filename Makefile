all:
	emcc wasm.cpp FFTOcean.cpp -o fftocean_wasm.js -s WASM=1 -O3 \
		-sMODULARIZE -sEXPORT_NAME=createModule \
		--bind -s ENVIRONMENT=web \
		-s EXPORT_ES6=1 \
		-s ALLOW_MEMORY_GROWTH=1 \
		-s EXPORTED_RUNTIME_METHODS=['HEAPF32'] \
		-msimd128