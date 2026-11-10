import sys
import importlib
import importlib.metadata as metadata

packages = [
    'pyscenic',
    'arboreto',
    'ctxcore',
    'pandas',
    'numpy',
    'loompy',
    'dask',
    'numba',
    'scanpy',
    'pyarrow'
]

print('python_executable\t' + str(sys.executable))
print('python_version\t' + str(sys.version).replace('\n', ' '))

for package in packages:
    try:
        module = importlib.import_module(package)
        try:
            version = metadata.version(package)
        except Exception:
            version = getattr(module, '__version__', 'installed_version_unknown')
        print(package + '\tPASS\t' + str(version) + '\t')
    except Exception as exc:
        text = str(exc).replace('\t', ' ').replace('\n', ' ')
        print(package + '\tFAIL\tNA\t' + text)
