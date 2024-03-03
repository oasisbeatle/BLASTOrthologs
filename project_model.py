from sqlobject import *
import os.path, sys

dbfile = 'proteinalignments.db3'

def init(new=False):
    # Magic formatting for database URI
    conn_str = os.path.abspath(dbfile)
    conn_str = 'sqlite:'+ conn_str
    # Connect to database
    sqlhub.processConnection = connectionForURI(conn_str)
    if new:
        # Create new tables (remove old ones if they exist)
        Proteins.dropTable(ifExists=True)
        Alignments.dropTable(ifExists=True)
        Proteins.createTable()
        Alignments.createTable()

class Proteins(SQLObject):
    accession = StringCol(alternateID=True)
    protein_name = StringCol()
    organism = StringCol()

class Alignments(SQLObject):
    query_protein = ForeignKey("Proteins")
    ref_protein = ForeignKey("Proteins")
    alignment_length = IntCol()
    alignment_score = FloatCol()
    eval = FloatCol()
    identity = IntCol()

