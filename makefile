PHONY: export_reqs

export_reqs:
	conda list -e > base_requirements.txt
	conda env export --no-builds > base_environment.yml

update_reqs:
	mamba activate base
	mamba env update --file base_environment.yml --prune

git_add:
	git status 
	git add -u 
	
install_reqs_wsl:
	sudo apt-get update
	sudo apt install build-essential