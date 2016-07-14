test:
	python3 tests/test_all.py

package:
	python3 setup.py bdist_wheel

readme_rst:
	grep -v "^\[!" README.md | sed -e "1d" > README.md.tmp
	pandoc --from=markdown --to=rst README.md.tmp -o README.rst
	rm README.md.tmp

md_to_html:
	pandoc README.md > README.html

clean:
	rm -rf README.html *~ build COVERnant.egg-info dist

pylint:
	pylint bin/bin covernantlib/* tests/*

pypi_submission:
	python3 setup.py bdist_wheel upload
	@echo "Go to https://pypi.python.org/pypi/COVERnant/"

new_release:
	@echo "* Create/checkout a release branch"
	@echo "  git branch release_v0.X"
	@echo "  git checkout release_v0.X"
	@echo "* Change bin/covernant"
	@echo "* Change setup.py"
	@echo "* Change CHANGELOG.txt"
	@echo "* Test package creation"
	@echo "* make package_to_pypi"
	@echo "* git add CHANGELOG.txt bin/covernant setup.py"
	@echo "* Commit changes e.g. 'git commit -m \"Set version to 0.X\"'"
	@echo "* Tag the commit e.g. 'git tag -a v0.X -m \"version v0.X\"'"
	@echo "* Merge release into dev and master"
	@echo "* Push it to github: git push"
	@echo "* Generate a new release based on this tag at"
	@echo "  https://github.com/konrad/COVERnant/releases/new"
