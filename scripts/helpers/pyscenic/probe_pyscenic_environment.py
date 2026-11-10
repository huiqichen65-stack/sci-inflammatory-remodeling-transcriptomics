import sys
import importlib
import importlib.metadata as metadata

packages = [
    'pyscenic',
    'arboreto',
    'ctxcore',
    'numpy',
    'pandas',
    'scipy',
    'sklearn',
    'numba',
    'dask',
    'distributed',
    'pyarrow',
    'loompy',
    'networkx'
]

print('PYTHON_EXECUTABLE\t' + sys.executable)
print('PYTHON_VERSION\t' + sys.version.replace('\n', ' '))

for package in packages:
    try:
        module = importlib.import_module(package)
        distribution_name = 'scikit-learn' if package == 'sklearn' else package
        try:
            version = metadata.version(distribution_name)
        except Exception:
            version = getattr(module, '__version__', 'UNKNOWN')
        print(package + '\tPASS\t' + str(version) + '\t')
    except Exception as exc:
        error = str(exc).replace('\t', ' ').replace('\n', ' ')
        print(package + '\tFAIL\tNA\t' + error)
