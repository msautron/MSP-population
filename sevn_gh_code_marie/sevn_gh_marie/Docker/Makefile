RKEY := $(shell git rev-parse --short HEAD)

# SEVNPY
sevnpy: sevnpy-amd64 sevnpy-arm64

sevnpy-push: sevnpy-amd64-push sevnpy-arm64-push

sevnpy-amd64: sevnpy.Dockerfile
	sudo docker build  --platform linux/amd64 --build-arg SEVN_VER=$(RKEY)  -t sevnpydocker:amd64 -f sevnpy.Dockerfile .

sevnpy-arm64: sevnpy.Dockerfile
	sudo docker build  --platform linux/arm64 --build-arg SEVN_VER=$(RKEY)  -t sevnpydocker:arm64 -f sevnpy.Dockerfile .

sevnpy-amd64-push: sevnpy-amd64
	sudo docker tag sevnpydocker:amd64 iogiul/sevnpydocker:amd64
	sudo docker push iogiul/sevnpydocker:amd64

sevnpy-arm64-push: sevnpy-arm64
	sudo docker tag sevnpydocker:arm64 iogiul/sevnpydocker:arm64
	sudo docker push iogiul/sevnpydocker:arm64

# SEVN
sevn: sevn-amd64 sevn-arm64

sevn-push: sevn-amd64-push sevn-arm64-push

sevn-amd64: sevn.Dockerfile
	sudo docker build  --platform linux/amd64 --build-arg SEVN_VER=$(RKEY)  -t sevndocker:amd64 -f sevn.Dockerfile .

sevn-arm64: sevn.Dockerfile
	sudo docker build  --platform linux/arm64 --build-arg SEVN_VER=$(RKEY)  -t sevndocker:arm64 -f sevn.Dockerfile .

sevn-amd64-push: sevn-amd64
	sudo docker tag sevndocker:amd64 iogiul/sevndocker:amd64
	sudo docker push iogiul/sevndocker:amd64

sevn-arm64-push: sevn-arm64
	sudo docker tag sevndocker:arm64 iogiul/sevndocker:arm64
	sudo docker push iogiul/sevndocker:arm64
