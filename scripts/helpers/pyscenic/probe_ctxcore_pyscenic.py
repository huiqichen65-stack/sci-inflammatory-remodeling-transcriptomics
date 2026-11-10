import sys
import importlib
import importlib.metadata as metadata

packages = ['setuptools', 'pkg_resources', 'ctxcore', 'pyscenic', 'arboreto', 'loompy']

print('python_executable\t' + sys.executable)
print('python_version\t' + sys.version.replace('\n', ' '))

for package in packages:
    try:
        module = importlib.import_module(package)
        distribution_name = 'setuptools' if package == 'pkg_resources' else package
        try:
            version = metadata.version(distribution_name)
        except Exception:
            version = getattr(module, '__version__', 'UNKNOWN')
        print(package + '\tPASS\t' + str(version) + '\t')
    except Exception as exc:
        text = str(exc).replace('\t', ' ').replace('\n', ' ')
        print(package + '\tFAIL\tNA\t' + text)
