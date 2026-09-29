"""Configuration failure tests, without mutating the real source submodule."""
import argparse
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cmake', default='cmake')
    parser.add_argument('--cc', required=True)
    parser.add_argument('--cxx', required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    with tempfile.TemporaryDirectory(prefix='STAR configure ') as temporary:
        root = Path(temporary)
        fixture = root / 'repository'
        for name in ['cmake', 'include/utility', 'src/utility', 'tests/star']:
            shutil.copytree(repo / name, fixture / name)

        def configure(name, source, options, diagnostic=None):
            result = subprocess.run([
                args.cmake, '-S', str(source), '-B', str(root / name),
                '-DCMAKE_C_COMPILER=' + str(Path(shutil.which(args.cc)).resolve()),
                '-DCMAKE_CXX_COMPILER=' + str(Path(shutil.which(args.cxx)).resolve()),
                '-DCMAKE_C_FLAGS=', '-DCMAKE_CXX_FLAGS=', '-DCMAKE_EXE_LINKER_FLAGS=',
                *options], text=True, capture_output=True, timeout=90)
            output = result.stdout + result.stderr
            if diagnostic:
                assert result.returncode != 0 and diagnostic in output, output
            else:
                assert result.returncode == 0, output

        configure('missing', fixture / 'tests/star', [], 'STAR source missing')
        configure('disabled', fixture / 'tests/star', ['-DRNANUE_BUILD_STAR=OFF'])
        fake = fixture / 'submodules/STAR'
        (fake / 'source').mkdir(parents=True)
        (fake / 'source/VERSION').write_text('#define STAR_VERSION "2.7.11b"\n')
        for command in [
            ['git', 'init', str(fake)], ['git', '-C', str(fake), 'add', '.'],
            ['git', '-C', str(fake), '-c', 'user.name=Test', '-c', 'user.email=test@example.invalid',
             '-c', 'commit.gpgsign=false', 'commit', '-m', 'deliberately wrong pin']]:
            subprocess.run(command, check=True, capture_output=True)
        configure('wrong-pin', fixture / 'tests/star', [], 'STAR must be pinned to')
        configure('missing-make', repo / 'tests/star',
                  ['-DRNANUE_STAR_MAKE=/nonexistent/gnu-make'], 'requires GNU Make')
        configure('missing-xxd', repo / 'tests/star',
                  ['-DRNANUE_STAR_XXD=/nonexistent/xxd'], 'requires a working xxd')
        configure('missing-openmp', repo / 'tests/star',
                  ['-DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE'], 'OpenMP')
        if platform.system() == 'Darwin':
            configure('static-mac', repo / 'tests/star', ['-DRNANUE_STAR_STATIC=ON'],
                      'RNANUE_STAR_STATIC requires Linux')
    print('Missing source, wrong pin, missing prerequisites and disabled configuration passed')


if __name__ == '__main__':
    main()
