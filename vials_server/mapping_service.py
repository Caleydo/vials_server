import os
import phovea_server

__author__ = 'Hendrik Strobelt'

config = phovea_server.config.view('vials_server')


class EnsemblGeneNameMapping():
  def __init__(self):
    self.from_id = 'ensembl'
    self.to_id = 'genename'
    self.mapping_db = os.path.join(config.mapping_dir, 'gene_names.sqlite')

  def get_sql_fields(self):
    """
    a method to return all parameter needed to access the sql fields

    :return: dbFile, dbName, fromField, toField, description
    """
    return self.mapping_db, 'GeneNames', 'EnsembleID', 'symbol', 'name'


class MappingService:
  def __init__(self):
    pass

  @staticmethod
  def get_handler(from_id, to_id):
    if from_id == 'ensembl' and to_id == 'genename':
      return EnsemblGeneNameMapping()

  @staticmethod
  def guess_id_type(samples):
    ensg_count = 0
    for sample in samples:
      if sample.startswith("ENSG"):
        ensg_count += 1

    if ensg_count > .8 * len(samples):
      return 'ensembl'
    else:
      return 'genename'
