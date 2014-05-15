DROP FUNCTION IF EXISTS FIND_COG;
CREATE FUNCTION FIND_COG(WID BIGINT(20)) RETURNS TEXT
RETURN (
SELECT IFNULL(GROUP_CONCAT(DISTINCT Id ORDER BY Id SEPARATOR ','),'-') FROM (
(SELECT g.Id AS Id
  FROM GeneBankCDS c INNER JOIN
  COGMember m ON c.Locus_Tag = m.Id INNER JOIN 
  COGOrthologousGroup g ON g.WID = m.COGOrthologousGroup_WID
  WHERE c.WID = WID)
UNION (
 SELECT COGId AS Id FROM GeneBankCOG WHERE GeneBankCDS_WID = WID
)) AS A
);

DROP FUNCTION IF EXISTS FIND_PROTCLUST;
CREATE FUNCTION FIND_PROTCLUST(WID BIGINT(20)) RETURNS TEXT
RETURN (SELECT IFNULL(GROUP_CONCAT(DISTINCT p.Entry ORDER BY p.Entry SEPARATOR ','),'-')
  FROM GeneBankCDS_has_GeneInfo cg INNER JOIN 
  ProtClust_has_GeneInfo pg ON pg.GeneInfo_WID = cg.GeneInfo_WID INNER JOIN
  ProtClust p ON p.WID = pg.ProtClust_WID WHERE GeneBankCDS_WID = WID);

SELECT 
    gb.Gi,
    gb.LocusName,
    gb.TaxId,
    IFNULL(c.ProteinId, '-'),
    IFNULL(c.Locus_Tag, '-'),
    l.pFrom,
    l.pTo,
    FIND_COG(c.WID),
    FIND_PROTCLUST(c.WID)
    INTO OUTFILE '/tmp/taxoner_DB.txt'
	LINES TERMINATED BY '\n'
FROM
    GeneBankCDS c
        INNER JOIN
    GeneBank gb ON c.GeneBank_WID = gb.WID
        INNER JOIN
    GeneBankCDSLocation l ON l.GeneBankCDS_WID = c.WID;


DROP FUNCTION IF EXISTS FIND_COG;
DROP FUNCTION IF EXISTS FIND_PROTCLUST;
