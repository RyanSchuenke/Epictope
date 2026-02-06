import subprocess
import os
import gzip
import shutil
import ftplib
import tempfile
from epictope.find_executable import find_exe

def ftp_download(species:str, cds_folder:os.PathLike, release: str = "release-108") -> os.PathLike:
    """
    Download the protein sequence files of a species from the uniprot ftp server
    
    :param species: species to download sequence data for
    :type species: str
    :param cds_folder: path to download sequences into
    :type cds_folder: os.PathLike
    :param release: version of uniprot to look for species data files
    :type release: str
    :return: path to downloaded species protein sequences file
    :rtype: PathLike
    """
    ftp = ftplib.FTP("ftp.ensembl.org")
    ftp.login()
    try:
        ftp.cwd("/pub/"+release+"/fasta/"+species+"/pep")
        for download_link in ftp.nlst():
            if download_link.endswith(".all.fa.gz"):
                ftp.retrbinary("RETR "+download_link, open(os.path.join(cds_folder,download_link), 'wb').write)
                break
    except Exception:
        ftp.quit()
        raise Exception("Could not find ftp link for species "+species)
    ftp.quit()
    return download_link

def make_protein_db(gz_file:str) -> None:
    """
    creates a protein database with the provided gzip file
    
    :param gz_file: file to make a database with
    :type gz_file: str
    """
    file = os.path.splitext(gz_file)[0]
    with gzip.open(gz_file, 'rb') as f_in:
        with open(file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    make_db_exe = find_exe("makeblastdb")
    subprocess.run([make_db_exe, "-in", file, "-dbtype", "prot", "-parse_seqids"])

def fetch_db(db:str, cds_folder:os.PathLike) -> os.PathLike:
    """
    Get the path to a database starting with the provided string
    
    :param db: name of database
    :type db: str
    :param cds_folder: folder to search for database in
    :type cds_folder: os.PathLike
    :return: path to database
    :rtype: os.PathLike
    """
    if not db:
        raise Exception("No BLAST database specified")
    for file in os.listdir(cds_folder):
        if os.path.isfile(os.path.join(cds_folder, file)) and file.lower().startswith(db) and file.endswith(".all.fa"):
            return os.path.join(cds_folder, file)
    else:
        raise Exception("Blast database for species "+db+" not found")

def install_db(species:list[str], cds_folder:os.PathLike, force:bool = False) -> None:
    """
    Make a blast database for all provided species
    
    :param species: list of species to make databases for
    :type species: list[str]
    :param cds_folder: folder to make databases in
    :type cds_folder: os.PathLike
    :param force: boolean to force remaking of databases if one already exists for the species
    :type force: bool
    """
    for s in species:
        try:
            if not force:
                fetch_db(s.lower(), cds_folder)
                print(f"Database for '{s}' already exists, skipping")
                continue
        except Exception:
            pass
        file = ftp_download(s.lower(), cds_folder)
        make_protein_db(os.path.join(cds_folder, file))
    print("all databases installed")

def fetch_seq(seq_id:str, db:str, cds_folder:os.PathLike, outfmt:str = "%s") -> str:
    """
    Retrieve a sequence from its database using the sequence id
    
    :param seq_id: sequence id to search blast database with
    :type seq_id: str
    :param db: database to search for sequence in
    :type db: str
    :param cds_folder: folder to find blast databases
    :type cds_folder: os.PathLike
    :param outfmt: outformat of sequence
    :type outfmt: str
    :return: sequence of the protein found in the database
    :rtype: str
    """
    dbcmd = find_exe("blastdbcmd")
    db_path = fetch_db(db, cds_folder)
    stdout = subprocess.check_output([dbcmd, '-db', db_path, '-entry', seq_id, '-outfmt', outfmt]).decode('utf-8')
    return stdout

def blast(seq:str, db: str, folders:list[os.PathLike], blast_type: str = "blastp", outfmt: str = '10 ', 
          custom_fmt:str = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore') -> list[list]:
    """
    Blast a sequence against a single database and return the results
    
    :param seq: sequence to use as query in blast
    :type seq: str
    :param db: database to blast query against
    :type db: str
    :param folders: dictionary of folders to find database and save output to
    :type folders: list[os.PathLike]
    :param blast_type: sequence type of blast being done
    :type blast_type: str
    :param outfmt: outformat specifier
    :type outfmt: str
    :param custom_fmt: columns to be included in output
    :type custom_fmt: str
    :return: list of blast hits
    :rtype: list[list]
    """
    db_path = fetch_db(db, folders["cds_folder"])
    
    blast_exe = find_exe(blast_type)
    
    outfmt += custom_fmt
    
    with tempfile.NamedTemporaryFile(mode='w', suffix=".fasta") as tmp:
        tmp.write(">query\n")
        tmp.write(seq+"\n")
        tmp.flush()
        
        stdout = subprocess.check_output([blast_exe, '-db', db_path, '-query', tmp.name, '-outfmt', outfmt]).decode('utf-8')
    
    out = []
    for l in stdout.splitlines():
        out.append(l.split(sep=','))
    
    return out
