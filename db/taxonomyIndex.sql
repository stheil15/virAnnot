PRAGMA encoding = "UTF-8";
PRAGMA default_synchronous = OFF;

CREATE INDEX tax_id ON names(tax_id);
CREATE INDEX unique_name ON names(unique_name);
CREATE INDEX name_txt ON names(name_txt);
CREATE INDEX parent_tax_id ON nodes(parent_tax_id);
CREATE INDEX prot_acc ON prot_accession2taxid(accession);
CREATE INDEX nucl_acc ON nucl_accession2taxid(accession);
