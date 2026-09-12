# A miniature population SV catalogue in the shape of the 1000 Genomes
# phase 3 integrated map, written to a tempdir so no catalogue file has to
# live in the repository.
mini_mei_vcf <- function(dir = tempdir(), gz = FALSE) {
  path <- fs::path(dir, if (gz) "mini_mei.vcf.gz" else "mini_mei.vcf")
  lines <- c(
    "##fileformat=VCFv4.1",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001",
    "6\t101194045\tSVA_umary_SVA_298\tT\t<INS:ME:SVA>\t.\t.\tTSD=AAAATAAAGAAACT;SVTYPE=SVA;SVLEN=267;AC=295;AF=0.05890575;AN=5008\tGT\t0|1",
    "6\t101205351\tALU_umary_ALU_5419\tA\t<INS:ME:ALU>\t.\t.\tTSD=null;SVTYPE=ALU;SVLEN=258;AC=56;AF=0.01118211;AN=5008\tGT\t0|0",
    "chr1\t645710\tL1_umary_LINE1_1\tA\t<INS:ME:LINE1>\t.\t.\tTSD=null;SVTYPE=LINE1;SVLEN=1200;AC=3;AF=0.0006;AN=5008\tGT\t0|0",
    "6\t51847748\tDEL_pindel_1\tA\t<CN0>\t.\tPASS\tSVTYPE=DEL;END=51848000;AC=10;AF=0.002;AN=5008\tGT\t0|0"
  )
  if (gz) {
    con <- gzfile(path, "w"); writeLines(lines, con); close(con)
  } else {
    write_lines(lines, path)
  }
  path
}

test_that("read_mei_catalogue parses INFO and keeps insertion classes only", {
  m <- read_mei_catalogue(mini_mei_vcf())
  expect_equal(nrow(m), 3)
  expect_setequal(m$svtype, c("SVA", "ALU", "LINE1"))
  expect_false("DEL" %in% m$svtype)

  sva <- m |> filter(svtype == "SVA")
  expect_equal(sva$pos, 101194045L)
  expect_equal(sva$tsd, "AAAATAAAGAAACT")
  expect_equal(sva$af, 0.05890575)
  expect_equal(sva$ac, 295L)
  expect_equal(sva$svlen, 267L)
  # An insertion record has no span, so end falls back to pos.
  expect_equal(sva$end, sva$pos)

  # TSD=null is a missing value, not the string "null".
  expect_true(is.na(m$tsd[m$svtype == "ALU"]))
  # Chromosomes come back unprefixed, as everywhere else in the toolkit.
  expect_equal(m$chrom[m$svtype == "LINE1"], "1")
})

test_that("read_mei_catalogue takes the class from the ALT when it is there", {
  # gnomAD puts every insertion under SVTYPE=INS and names the element class
  # in the ALT allele, so reading SVTYPE alone would collapse all three into
  # one class and lose the comparison the survey turns on.
  path <- fs::path(tempdir(), "mini_gnomad.vcf")
  write_lines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "1\t100\tgnomAD-SV_v2.1_ALU_1\tN\t<INS:ME:ALU>\t100\tPASS\tEND=101;SVTYPE=INS;SVLEN=300;AC=145;AN=21366;AF=0.006786",
    "1\t200\tgnomAD-SV_v2.1_SVA_1\tN\t<INS:ME:SVA>\t100\tPASS\tEND=201;SVTYPE=INS;SVLEN=1200;AC=10;AN=21366;AF=0.00047",
    "1\t300\tgnomAD-SV_v2.1_INS_1\tN\t<INS>\t100\tPASS\tEND=301;SVTYPE=INS;SVLEN=50;AC=5;AN=21366;AF=0.0002"), path)
  m <- read_mei_catalogue(path)
  expect_equal(nrow(m), 2)
  expect_setequal(m$svtype, c("ALU", "SVA"))
  # The untyped <INS> record is not a mobile element and is dropped.
  expect_false(any(str_detect(m$id, "INS_1")))
  expect_equal(m$af[m$svtype == "ALU"], 0.006786)
})

test_that("read_mei_catalogue reads gzip and honours types = NULL", {
  expect_equal(nrow(read_mei_catalogue(mini_mei_vcf(gz = TRUE))), 3)
  all_rec <- read_mei_catalogue(mini_mei_vcf(), types = NULL)
  expect_equal(nrow(all_rec), 4)
  # The deletion record carries an END, which is how its size is recovered.
  del <- all_rec |> filter(svtype == "DEL")
  expect_equal(del$end - del$pos, 252L)
  expect_equal(nrow(read_mei_catalogue(mini_mei_vcf(), types = "SVA")), 1)
  expect_error(read_mei_catalogue(fs::path(tempdir(), "nope.vcf")), "not found")
})

test_that("breakends_near_mei tests both ends and respects the window", {
  mei <- read_mei_catalogue(mini_mei_vcf())
  sv <- tibble(
    UUID = c("u1", "u2", "u3"),
    Sample = c("S1", "S1", "S2"),
    sv_class = c("DEL", "TRA", "DEL"),
    CHROM_A = c("6", "1", "6"),
    START_A = c(51847748L, 645760L, 90000000L),
    CHROM_B = c("6", "6", "6"),
    START_B = c(101194032L, 101205351L, 95000000L)
  )
  h <- breakends_near_mei(sv, mei, window = 500L)
  # u1 matches at its B end only: the A end is on a deletion record, which
  # is not an insertion and was filtered out of the catalogue.
  expect_equal(sort(unique(h$UUID)), c("u1", "u2"))
  expect_equal(h$end[h$UUID == "u1"], "B")
  expect_equal(h$mei_id[h$UUID == "u1"], "SVA_umary_SVA_298")
  expect_equal(h$distance_bp[h$UUID == "u1"], 13L)
  # u2 matches at both ends, on two different chromosomes.
  expect_setequal(h$end[h$UUID == "u2"], c("A", "B"))
  # u3 is 5 Mb from anything.
  expect_false("u3" %in% h$UUID)

  # A tighter window drops the 13 bp offset only when it is tighter than 13.
  expect_true("u1" %in% breakends_near_mei(sv, mei, window = 13L)$UUID)
  expect_false("u1" %in% breakends_near_mei(sv, mei, window = 12L)$UUID)

  expect_error(breakends_near_mei(sv |> select(-START_B), mei), "missing columns")
  expect_error(breakends_near_mei(sv, mei |> select(-svtype)), "missing columns")
})

test_that("breakends_near_mei is not confused by the catalogue's end column", {
  # The catalogue carries `end` and the result labels breakends `end` too.
  # Joining without dropping one makes the bare symbol resolve to base::end.
  mei <- read_mei_catalogue(mini_mei_vcf())
  expect_true("end" %in% names(mei))
  sv <- tibble(UUID = "u1", Sample = "S1", sv_class = "DEL", CHROM_A = "6",
               START_A = 101194045L, CHROM_B = "6", START_B = 1L)
  h <- breakends_near_mei(sv, mei)
  expect_equal(h$end, "A")
  expect_type(h$end, "character")
})

test_that("breakends_near_mei survives a cohort with no sv_class column", {
  mei <- read_mei_catalogue(mini_mei_vcf())
  sv <- tibble(UUID = "u1", Sample = "S1", CHROM_A = "6", START_A = 101194045L,
               CHROM_B = "6", START_B = 1L)
  h <- breakends_near_mei(sv, mei)
  expect_equal(nrow(h), 1)
  expect_true(is.na(h$sv_class))
})

test_that("positions_in_filter reads BED and both ends of a BEDPE", {
  bed <- fs::path(tempdir(), "f.bed")
  write_lines(c("6\t101194000\t101195000", "1\t1000\t2000"), bed)
  bedpe <- fs::path(tempdir(), "f.bedpe")
  write_lines("6\t51847000\t51848000\t6\t900000\t901000\tpair\t0\t+\t-", bedpe)

  pos <- tibble(chrom = c("6", "6", "6", "2"),
                pos = c(101194045L, 51847748L, 900500L, 5L))
  expect_equal(positions_in_filter(pos, bed), c(TRUE, FALSE, FALSE, FALSE))
  # Both ends of the BEDPE are treated as independent intervals.
  expect_equal(positions_in_filter(pos, bedpe), c(FALSE, TRUE, TRUE, FALSE))
  # slop widens the intervals.
  expect_true(positions_in_filter(tibble(chrom = "1", pos = 2500L), bed, slop = 1000L))
  expect_false(positions_in_filter(tibble(chrom = "1", pos = 2500L), bed))

  expect_error(positions_in_filter(pos, fs::path(tempdir(), "nope.bed")), "not found")
  expect_error(positions_in_filter(tibble(x = 1), bed), "missing columns")
})
