import pathlib
import subprocess

def download_git_repo(version, url, dest_dir):
    """Function to download a git repo"""
    # Delete old output dir if existing and create parent dirs if not existing
    try:
        rm_dir = subprocess.run(['rm','-rf', str(dest_dir)], check = True, timeout = 60)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as err:
        rm_dir.kill()
        raise

    if not isinstance(dest_dir, pathlib.PosixPath):
        dest_dir = pathlib.Path(dest_dir)

    dest_dir.parent.mkdir(exist_ok = True)
    # Download
    try:
        downloading = subprocess.run(['git', 'clone', 
                                        '-b', version, 
                                        '--single-branch', '--depth=1', 
                                        url, str(dest_dir)],
                                        check = True,
                                        timeout = 500)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as err:
        downloading.kill()
        raise
        

def get_commit_git(gitrepo_dir):
    """Function to get the commit number from a folder (must be a git repo)"""
    try:
        commit = subprocess.check_output(['git', 
                                        '--git-dir', 
                                        '{}/.git'.format(str(gitrepo_dir)), 
                                        'log', '-n', '1', '--pretty=format:"%H"'],
                                        timeout = 30)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as err:
        commit.kill()
        raise                    
    return commit

def download_software_kmerfinder(kmerfinder_software_dir, version = '3.0.2'):
    """Function to download kmerfinder if it is not present"""
    if not isinstance(kmerfinder_software_dir, pathlib.PosixPath):
        kmerfinder_software_dir = pathlib.Path(kmerfinder_software_dir)
    if not kmerfinder_software_dir.joinpath('kmerfinder.py').is_file():
        print("\x1b[0;33m Downloading kmerfinder software...\n\033[0;0m")
        download_git_repo(version, 
                        'https://bitbucket.org/genomicepidemiology/kmerfinder.git',
                        kmerfinder_software_dir)
    return version
    
def download_software_mlst7(mlst7_software_dir, version = '2.0.4'):
    """Function to download MLST (CGE) if it is not present"""
    if not isinstance(mlst7_software_dir, pathlib.PosixPath):
        mlst7_software_dir = pathlib.Path(mlst7_software_dir)
    if not mlst7_software_dir.joinpath('mlst.py').is_file():
        print("\x1b[0;33m Downloading MLST (CGE) software...\n\033[0;0m")
        download_git_repo(version, 
                        'https://bitbucket.org/genomicepidemiology/mlst.git',
                        mlst7_software_dir)
    return version

def download_db_kmerfinder(kmerfinder_db_dir, version = '20210425'):
    """Function to download kmerfinder database if it is not present"""
    if not isinstance(kmerfinder_db_dir, pathlib.PosixPath):
        kmerfinder_db_dir = pathlib.Path(kmerfinder_db_dir)
    if not kmerfinder_db_dir.joinpath('config').exists():
        print("\x1b[0;33m Downloading KmerFinder database...\n\033[0;0m")
        download_git_repo('master', 
                        'https://bitbucket.org/genomicepidemiology/kmerfinder_db.git',
                        kmerfinder_db_dir)
        try:
            build = subprocess.run(['bash', 
                                'INSTALL.sh',
                                str(kmerfinder_db_dir.absolute()), 
                                'bacteria', 
                                version],
                                cwd = str(kmerfinder_db_dir),
                                check=True,
                                timeout = 3000)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as err:
            build.kill()
            raise
        version = get_commit_git(kmerfinder_db_dir)
    return version

def download_db_mlst7(mlst7_db_dir, version = 'master'):
    """Function to download the MLST (CGE) database if it is not present"""
    if not isinstance(mlst7_db_dir, pathlib.PosixPath):
        mlst7_db_dir = pathlib.Path(mlst7_db_dir)
    if not mlst7_db_dir.joinpath('senterica', 'senterica.length.b').is_file():
        print("\x1b[0;33m Downloading MLST (CGE) database...\n\033[0;0m")
        download_git_repo(version, 
                        'https://bitbucket.org/genomicepidemiology/mlst_db.git',
                        mlst7_db_dir)
        try:
            build = subprocess.run(['python', 
                            'INSTALL.py',
                            'kma_index'], 
                            check = True,
                            cwd = str(mlst7_db_dir),
                            timeout=800)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as err:
            build.kill()
            raise 
        version = get_commit_git(mlst7_db_dir)
    return version

def download_db_serotypefinder(serotypefinder_db_dir, version = 'master'):
    """Function to download the SerotypeFinder database if it is not present"""
    if not isinstance(serotypefinder_db_dir, pathlib.PosixPath):
        serotypefinder_db_dir = pathlib.Path(serotypefinder_db_dir)
    if not serotypefinder_db_dir.joinpath('H_type.seq.b').is_file():
        print("\x1b[0;33m Downloading SerotypeFinder database...\n\033[0;0m")
        download_git_repo(version, 
                        'https://bitbucket.org/genomicepidemiology/serotypefinder_db.git',
                        serotypefinder_db_dir)
        try:
            build = subprocess.run(['python', 
                            str('INSTALL.py'),
                            'kma_index'], 
                            cwd = str(serotypefinder_db_dir),
                            check=True,
                            timeout = 800)        
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as err:
            build.kill()
            raise 
        version = get_commit_git(serotypefinder_db_dir)
    return version


def download_db_seroba(seroba_db_dir, version='master', kmersize = 71):
    """Function to download the Seroba database if it is not present"""
    if not isinstance(seroba_db_dir, pathlib.PosixPath):
        seroba_db_dir = pathlib.Path(seroba_db_dir)
    if not seroba_db_dir.joinpath('database', 'cdhit_cluster').is_file():
        print("\x1b[0;33m Downloading Seroba database...\n\033[0;0m")
        download_git_repo(version, 
                        'https://github.com/sanger-pathogens/seroba.git',
                        seroba_db_dir)
        try:
            kmersize = str(kmersize)
            rm_dir = subprocess.run(['rm', '-rf', str(seroba_db_dir.joinpath('scripts')),
                                            '&&', 'rm', '-rf', str(seroba_db_dir.joinpath('seroba'))],
                                        check = True, timeout=60)
            build = subprocess.run(['seroba', 'createDBs', 
                                    'database',
                                    kmersize], 
                                    cwd = str(seroba_db_dir),
                                    check=True,
                                    timeout = 800)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as err:
            build.kill()
            raise 
        version = get_commit_git(seroba_db_dir)
    return version

def get_downloads_juno_typing(db_dir, current_dir, update_dbs, seroba_kmersize):
    if not isinstance(db_dir, pathlib.PosixPath):
        db_dir = pathlib.Path(db_dir)
    if not isinstance(current_dir, pathlib.PosixPath):
        current_dir = pathlib.Path(current_dir)
    if update_dbs:
        try:
            rm_dir = subprocess.run(['rm', '-rf', str(db_dir)], check = True, timeout = 60)
        except:
            rm_dir.kill()
            raise
    software_version = {'kmerfinder': download_software_kmerfinder(current_dir.joinpath('kmerfinder')),
                        'mlst7': download_software_mlst7(current_dir.joinpath('cge-mlst')),
                        'kmerfinder_db': download_db_kmerfinder(db_dir.joinpath('kmerfinder_db')),
                        'mlst7_db': download_db_mlst7(db_dir.joinpath('mlst7_db')),
                        'serotypefinder_db': download_db_serotypefinder(db_dir.joinpath('serotypefinder_db')),
                        'seroba_db': download_db_seroba(db_dir.joinpath('seroba_db'), kmersize = seroba_kmersize)}
    return software_version

# if __name__ == '__main__':
#     current_dir = pathlib.Path(__file__).parent.absolute()
#     db_dir = current_dir.joinpath('fake_db')
#     update_dbs = False
#     seroba_kmersize=71
#     versions = get_downloads_juno_typing(db_dir, current_dir, update_dbs, seroba_kmersize)
#     print(versions)