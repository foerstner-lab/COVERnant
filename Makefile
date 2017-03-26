test:
	python -m pytest

package:
	python3 setup.py bdist

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
	python3 setup.py bdist upload
	@echo "Go to https://pypi.python.org/pypi/COVERnant/"

new_release:
	@echo "* Create/checkout a release branch"
	@echo "  $$ git branch release_v0.X.0"
	@echo "  $$ git checkout release_v0.X.0"
	@echo "* Change bin/covernant"
	@echo "* Change setup.py"
	@echo "* Change CHANGELOG.txt"
	@echo "* Test package creation"
	@echo "  $$ make package"
	@echo "* Submit package"
	@echo "  $$ make package_to_pypi"
	@echo "* Add and commit changes"
	@echo "  $$ git add CHANGELOG.txt bin/covernant setup.py"
	@echo "  $$ git commit -m \"Set version to 0.X.0\""
	@echo "* Tag the commit"
	@echo "  $$ git tag -a v0.X.0 -m \"version v0.X.0\""
	@echo "* Merge release into dev and master and push"
	@echo " $ git checkout master"
	@echo " $ git merge release_v0.X.0"
	@echo " $ git push"
	@echo " $ git checkout dev"
	@echo " $ git merge release_v0.X.0"
	@echo " $ git push"
	@echo "* Generate a new release based on this tag at"
	@echo "  https://github.com/konrad/COVERnant/releases/new"

test_package:
	@echo "$$ make package"
	@echo "$$ virtualenv covernant_test"
	@echo "$$ cd covernant_test"
	@echo "$$ source bin/activate"
	@echo "$$ pip3 install ../dist/COVERnant-*-py3-none-any.whl"
	@echo "$$ bin/covernant"
