"""Check toolchain notice discovery without depending on the host's layout."""
import argparse
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cmake', required=True)
    args = parser.parse_args()
    module = Path(__file__).resolve().parents[2] / 'cmake/star_runtime_license.cmake'
    # Synthetic markers exercise discovery and validation, not licence wording.
    notice = 'GCC RUNTIME LIBRARY EXCEPTION\nfixture\n2. No Weakening of GCC Copyleft.\n'
    with tempfile.TemporaryDirectory(prefix='STAR runtime notice ') as temporary:
        root = Path(temporary)
        script = root / 'stage.cmake'
        script.write_text(f'include([[{module}]])\nrnanue_stage_star_runtime_license("${{destination}}")\n')

        def stage(name, relative=None, contents=notice, override=None, diagnostic=None):
            prefix = root / name / 'toolchain'
            compiler = prefix / 'bin/g++-14'
            compiler.parent.mkdir(parents=True)
            compiler.touch()
            if relative:
                source = prefix / relative
                source.parent.mkdir(parents=True, exist_ok=True)
                source.write_text(contents)
            output = root / name / 'build/notices/GCC-RUNTIME-EXCEPTION.txt'
            command = [args.cmake, f'-DCMAKE_CXX_COMPILER={compiler}',
                       '-DCMAKE_CXX_COMPILER_VERSION=14.3.0', f'-Ddestination={output}']
            if override is not None:
                command.append(f'-DRNANUE_STAR_GCC_RUNTIME_LICENSE={override}')
            result = subprocess.run(command + ['-P', str(script)], capture_output=True,
                                    text=True, timeout=15)
            text = result.stdout + result.stderr
            if diagnostic:
                assert result.returncode and diagnostic in text, text
                assert not output.exists()
            else:
                assert result.returncode == 0, text
                assert output.read_text() == contents
            # Staging must never create a second copy in the toolchain/source tree.
            assert not (prefix / 'GCC-RUNTIME-EXCEPTION.txt').exists()

        stage('homebrew layout', 'COPYING.RUNTIME')
        stage('debian layout', 'share/doc/gcc-14-base/copyright',
              'Distribution notices\n' + notice + 'Additional notices\n')
        stage('fedora layout', 'share/licenses/libgcc/COPYING.RUNTIME')
        custom = root / 'custom copyright'
        custom.write_text(notice)
        stage('explicit override', override=custom)
        stage('missing', diagnostic='Cannot find the GCC runtime exception')
        stage('missing override', override=root / 'absent',
              diagnostic='Cannot find the GCC runtime exception')
        custom.write_text('not a runtime exception')
        stage('invalid override', override=custom,
              diagnostic='does not contain the GCC runtime exception')
    print('Runtime notice discovery, byte-preserving staging and diagnostics passed')


if __name__ == '__main__':
    main()
