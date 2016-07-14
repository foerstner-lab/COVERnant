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
	rm -f README.html *~
