-- Provide statistics on what fraction of seed alignments are Cyanobacterial per family

SELECT
    *,
    (res.num_cyano / res.num_bact) AS proportion
FROM (
    -- number of bacterial and  cyanobacterial seeds per family
    SELECT
        sr.rfam_acc, fam.rfam_id, fam.description, fam.num_seed,
        COUNT(*) AS num_bact,
        SUM(is_cyano) AS num_cyano
    FROM seed_region as sr
    INNER JOIN (
        -- Indicate which bacterial rfamseq entries are cyanobacterial
        SELECT DISTINCT
            rs.rfamseq_acc,
            tax_string LIKE 'Bacteria; Cyanobacteria;%' AS is_cyano
        FROM taxonomy AS tax
        INNER JOIN rfamseq AS rs ON rs.ncbi_id = tax.ncbi_id
        -- subset to bacterial entries
        WHERE tax.tax_string LIKE 'Bacteria;%'
    ) AS  myrs ON myrs.rfamseq_acc = sr.rfamseq_acc
    INNER JOIN family AS fam ON fam.rfam_acc = sr.rfam_acc
    GROUP BY sr.rfam_acc, fam.rfam_id, fam.description, fam.num_seed
) AS res
ORDER BY proportion DESC;
