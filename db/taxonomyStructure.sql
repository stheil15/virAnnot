CREATE TABLE `prot_accession2taxid` (
  `accession` varchar(255) NOT NULL DEFAULT '',
  `accession.version` varchar(255) NOT NULL DEFAULT '',
  `taxid` double NOT NULL,
  `gi` double NOT NULL,
  PRIMARY KEY (`gi`)
);
CREATE TABLE `nucl_accession2taxid` (
  `accession` varchar(255) NOT NULL DEFAULT '',
  `accession.version` varchar(255) NOT NULL DEFAULT '',
  `taxid` double NOT NULL,
  `gi` double NOT NULL,
  PRIMARY KEY (`gi`)
);
CREATE TABLE `names` (
  `tax_id` int unsigned NOT NULL DEFAULT '0',
  `name_txt` varchar(255) NOT NULL DEFAULT '',
  `unique_name` varchar(255) DEFAULT NULL,
  `name_class` varchar(32) NOT NULL DEFAULT ''
);
CREATE TABLE `nodes` (
  `tax_id` int unsigned NOT NULL DEFAULT '0',
  `parent_tax_id` int unsigned NOT NULL DEFAULT '0',
  `rank` varchar(32) DEFAULT NULL,
  `embl_code` varchar(16) DEFAULT NULL,
  `division_id` int NOT NULL DEFAULT '0',
  `inherited_div_flag` int NOT NULL DEFAULT '0',
  `genetic_code_id` int NOT NULL DEFAULT '0',
  `inherited_GC_flag` int NOT NULL DEFAULT '0',
  `mitochondrial_genetic_code_id` int NOT NULL DEFAULT '0',
  `inherited_MGC_flag` int NOT NULL DEFAULT '0',
  `GenBank_hidden_flag` int NOT NULL DEFAULT '0',
  `hidden_subtree_root_flag` int NOT NULL DEFAULT '0',
  `comments` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`tax_id`)
);
CREATE TABLE `cdd_taxo` (
  `id` int unsigned NOT NULL DEFAULT '0',
  `family` varchar(32) DEFAULT NULL,
  `genus` varchar(32) DEFAULT NULL,
  `no rank` varchar(32) DEFAULT NULL,
  `specie` varchar(32) DEFAULT NULL,
  `superkingdom` varchar(32) DEFAULT NULL
);
CREATE INDEX tax_id ON names(tax_id);
CREATE INDEX unique_name ON names(unique_name);
CREATE INDEX name_txt ON names(name_txt);
CREATE INDEX parent_tax_id ON nodes(parent_tax_id);
CREATE INDEX prot_acc ON prot_accession2taxid(accession);
CREATE INDEX nucl_acc ON nucl_accession2taxid(accession);
CREATE INDEX nucl_taxid ON nucl_accession2taxid(taxid);
CREATE INDEX prot_taxid ON nucl_accession2taxid(taxid);
CREATE TABLE `virus2hosts` (
`representative` varchar(32) DEFAULT NULL,
`neighbor` varchar(32) DEFAULT NULL,
`taxonomy` varchar(32) DEFAULT NULL,
`lineage` varchar(32) DEFAULT NULL,
`host` varchar(32) DEFAULT NULL
);
