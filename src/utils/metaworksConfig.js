const fs = require("fs");
const path = require("path");
const yaml = require("js-yaml");

const PROTEIN_MARKERS = [
  "COI",
  "rbcL_eukaryota",
  "rbcL_diatom",
  "rbcL_landPlant",
];

const PSEUDOGENE_DEFAULTS = {
  COI: {
    removal_type: "1",
    genetic_code: "5",
    grep_type: "1",
    taxon1: "Arthropoda",
  },
  rbcL_landPlant: {
    removal_type: "1",
    genetic_code: "1",
    grep_type: "1",
    taxon1: "Embryophyta",
  },
  rbcL_diatom: {
    removal_type: "1",
    genetic_code: "1",
    grep_type: "1",
    taxon1: "Bacillariophyta",
  },
  rbcL_eukaryota: {
    removal_type: "1",
    genetic_code: "1",
    grep_type: "1",
    taxon1: "",
  },
};

function applyPseudogeneDefaults(services, marker) {
  const defaults = PSEUDOGENE_DEFAULTS[marker];
  if (!defaults) return false;
  const list = Array.isArray(services) ? services : [services];
  const panel = list.find((service) => service.serviceName === "pseudogene filtering");
  if (!panel) return false;
  panel.selected = true;
  (panel.Inputs || []).forEach((field) => {
    if (Object.prototype.hasOwnProperty.call(defaults, field.name)) {
      field.value = defaults[field.name];
    }
  });
  return true;
}

const IUPAC_COMPLEMENT = {
  A: "T",
  T: "A",
  C: "G",
  G: "C",
  R: "Y",
  Y: "R",
  S: "S",
  W: "W",
  K: "M",
  M: "K",
  B: "V",
  V: "B",
  D: "H",
  H: "D",
  N: "N",
};

class MetaWorksConfigError extends Error {
  constructor(message, code) {
    super(message);
    this.name = "MetaWorksConfigError";
    this.code = code;
  }
}

function setting(services, name) {
  const list = Array.isArray(services) ? services : [services];
  for (let i = 0; i < list.length; i++) {
    const inputs = [].concat(list[i].Inputs || [], list[i].extraInputs || []);
    const found = inputs.find((input) => input.name === name);
    if (found) return found.value;
  }
  throw new MetaWorksConfigError(`Missing MetaWorks setting: ${name}`, "FORM");
}

function needsClassifier(marker) {
  return marker !== "16S" && marker !== "28S_fungi";
}

function reverseComplement(sequence) {
  const bases = String(sequence).trim().toUpperCase().split("");
  if (bases.length === 0) {
    throw new MetaWorksConfigError("A primer sequence is empty.", "FORM");
  }
  const complemented = bases.map((base) => {
    const partner = IUPAC_COMPLEMENT[base];
    if (!partner) {
      throw new MetaWorksConfigError(
        `Primer "${sequence}" contains "${base}", which is not an IUPAC DNA base.`,
        "FORM"
      );
    }
    return partner;
  });
  return complemented.reverse().join("");
}

function primerList(value, label) {
  const primers = (Array.isArray(value) ? value : [value])
    .map((primer) => String(primer).trim())
    .filter((primer) => primer.length > 0);
  if (primers.length === 0) {
    throw new MetaWorksConfigError(`Add at least one ${label} primer.`, "FORM");
  }
  return primers;
}

function buildAdaptersFasta(forwardPrimers, reversePrimers) {
  const forward = primerList(forwardPrimers, "forward");
  const reverse = primerList(reversePrimers, "reverse");
  if (forward.length !== reverse.length) {
    throw new MetaWorksConfigError(
      "Add the same number of forward and reverse primers. Each pair becomes one amplicon.",
      "FORM"
    );
  }
  return forward
    .map((fwd, index) => {
      const rev = reverseComplement(reverse[index]);
      reverseComplement(fwd);
      return `>amplicon${index + 1};\n^${fwd.toUpperCase()}...${rev}$`;
    })
    .join("\n") + "\n";
}

function listTopLevelFastqGz(dirPath) {
  let entries;
  try {
    entries = fs.readdirSync(dirPath, { withFileTypes: true });
  } catch (error) {
    throw new MetaWorksConfigError(
      `Could not read the selected folder: ${error.message}`,
      "FLAT"
    );
  }
  return entries
    .filter((entry) => entry.isFile())
    .map((entry) => entry.name)
    .filter((name) => /\.(fastq|fq)\.gz$/i.test(name))
    .sort();
}

function splitMate(filename) {
  const extMatch = filename.match(/^(.*)(\.(?:fastq|fq)\.gz)$/i);
  if (!extMatch) {
    return null;
  }
  const stem = extMatch[1];
  const ext = extMatch[2];
  const readMarker = /([._])R([12])(?![A-Za-z0-9])/gi;
  let found = null;
  let match;
  while ((match = readMarker.exec(stem)) !== null) {
    found = match;
  }
  if (!found) {
    return null;
  }
  return {
    head: stem.slice(0, found.index),
    sep: found[1] + found[0].charAt(1),
    read: found[2],
    tail: stem.slice(found.index + found[0].length),
    ext,
  };
}

function inferReadPattern(files) {
  if (files.length === 0) {
    throw new MetaWorksConfigError(
      "MetaWorks needs gzipped paired-end fastq files (.fastq.gz or .fq.gz) directly in the selected folder. Subfolders are not used.",
      "FLAT"
    );
  }

  const parsed = files.map((filename) => ({ filename, mate: splitMate(filename) }));
  const unmarked = parsed.filter((item) => !item.mate).map((item) => item.filename);
  if (unmarked.length > 0) {
    throw new MetaWorksConfigError(
      `Could not find an R1/R2 marker in: ${unmarked.slice(0, 5).join(", ")}. Enter a pattern such as {sample}_L001_R{read}_001.fastq.gz`,
      "PATTERN"
    );
  }

  const emptyHeads = parsed.filter((item) => item.mate.head.length === 0);
  if (emptyHeads.length > 0) {
    throw new MetaWorksConfigError(
      "A filename starts with the read marker, so the sample name is missing.",
      "PATTERN"
    );
  }

  const signatures = new Set(
    parsed.map((item) => `${item.mate.sep}|${item.mate.tail}|${item.mate.ext}`)
  );
  if (signatures.size !== 1) {
    throw new MetaWorksConfigError(
      "These files do not share one R1/R2 filename pattern. Enter a pattern such as {sample}_L001_R{read}_001.fastq.gz",
      "PATTERN"
    );
  }

  const readsBySample = new Map();
  parsed.forEach((item) => {
    if (!readsBySample.has(item.mate.head)) {
      readsBySample.set(item.mate.head, new Set());
    }
    readsBySample.get(item.mate.head).add(item.mate.read);
  });

  const incomplete = [];
  readsBySample.forEach((reads, head) => {
    if (!reads.has("1") || !reads.has("2") || reads.size !== 2) {
      incomplete.push(head);
    }
  });
  if (incomplete.length > 0) {
    throw new MetaWorksConfigError(
      `Each sample needs one R1 file and one R2 file. Check: ${incomplete.slice(0, 5).join(", ")}`,
      "FLAT"
    );
  }

  const { sep, tail, ext } = parsed[0].mate;
  return `{sample}${sep}{read}${tail}${ext}`;
}

function patternToRegExp(pattern) {
  const sampleToken = "\u0000SAMPLE\u0000";
  const readToken = "\u0000READ\u0000";
  const tokenized = pattern
    .split("{sample}")
    .join(sampleToken)
    .split("{read}")
    .join(readToken);
  const escaped = tokenized.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
  const source = escaped
    .split(sampleToken)
    .join("(.+)")
    .split(readToken)
    .join("([12])");
  return new RegExp(`^${source}$`);
}

function assertPatternMatches(files, pattern) {
  if (pattern.includes("/") || pattern.includes("\\")) {
    throw new MetaWorksConfigError(
      "The filename pattern is only the file name, not a directory.",
      "PATTERN"
    );
  }
  if (!pattern.includes("{sample}") || !pattern.includes("{read}")) {
    throw new MetaWorksConfigError(
      "The pattern must contain {sample} and {read}.",
      "PATTERN"
    );
  }
  const expression = patternToRegExp(pattern);
  const unmatched = files.filter((filename) => !expression.test(filename));
  if (unmatched.length > 0) {
    throw new MetaWorksConfigError(
      `That pattern does not match: ${unmatched.slice(0, 5).join(", ")}`,
      "PATTERN"
    );
  }
}

function asNumber(value, label) {
  const number = Number(value);
  if (!Number.isFinite(number)) {
    throw new MetaWorksConfigError(`${label} must be a number.`, "FORM");
  }
  return number;
}

function resourcesFromManager(memoryBytes, cpuCount) {
  const bytes = Number(memoryBytes);
  const cpus = Number(cpuCount);
  const memoryGb = Number.isFinite(bytes) && bytes > 0
    ? Math.max(1, Math.floor(bytes / (1024 ** 3)))
    : 1;
  const threads = Number.isFinite(cpus) && cpus >= 1 ? Math.floor(cpus) : 1;
  return { memoryGb, threads };
}

function taxonName(value) {
  return String(value ?? "").trim().replace(/^-[ev]\s+/i, "").trim();
}

function panelSelected(services, serviceName) {
  const list = Array.isArray(services) ? services : [services];
  const panel = list.find((service) => service.serviceName === serviceName);
  return Boolean(panel && panel.selected === true);
}

function buildConfigObject(services, readPattern, classifierContainerPath, resources) {
  const marker = setting(services, "marker");
  const protein = PROTEIN_MARKERS.includes(marker);
  const filterPseudogenes = protein && panelSelected(services, "pseudogene filtering");
  const removalType = Number(setting(services, "removal_type"));
  if (filterPseudogenes && removalType === 2 && marker !== "COI") {
    throw new MetaWorksConfigError(
      "HMM pseudogene filtering (removal type 2) is only available for COI. bold.hmm is an arthropod COI profile.",
      "FORM"
    );
  }
  const grepType = Number(setting(services, "grep_type"));
  const kept = taxonName(setting(services, "taxon1"));
  const dropped = taxonName(setting(services, "taxon2"));
  if (filterPseudogenes && kept === "") {
    throw new MetaWorksConfigError(
      "Type the group to keep, for example Chlorophyta.",
      "FORM"
    );
  }
  if (filterPseudogenes && grepType === 2 && dropped === "") {
    throw new MetaWorksConfigError(
      "Type the group to discard, for example Chordata.",
      "FORM"
    );
  }

  const { memoryGb, threads } = resources;

  const custom = needsClassifier(marker) ? "yes" : "no";
  const classifierPath = custom === "yes" ? classifierContainerPath : "/dev/null";

  return {
    raw: "/input",
    raw_sample_read_wildcards: `/input/${readPattern}`,
    raw_sample_forward_wildcard: `/input/${readPattern.replace("{read}", "1")}`,
    raw_sample_reverse_wildcard: `/input/${readPattern.replace("{read}", "2")}`,
    dir: "/input/metaworks_out",
    SEQPREP: {
      q: asNumber(setting(services, "seqprep_q"), "SeqPrep quality"),
      o: asNumber(setting(services, "seqprep_o"), "SeqPrep overlap"),
      m: asNumber(setting(services, "seqprep_m"), "SeqPrep mismatch fraction"),
      n: asNumber(setting(services, "seqprep_n"), "SeqPrep match fraction"),
    },
    CUTADAPT: {
      fasta: "/input/metaworks_out/adapters.fasta",
      m: asNumber(setting(services, "cutadapt_m"), "minimum length"),
      q: String(setting(services, "cutadapt_q")),
      e: asNumber(setting(services, "cutadapt_e"), "Cutadapt error rate"),
      O: asNumber(setting(services, "cutadapt_O"), "adapter overlap"),
      mn: asNumber(setting(services, "cutadapt_mn"), "maximum Ns"),
      rc: String(setting(services, "cutadapt_rc")),
    },
    VSEARCH_DENOISE: {
      minsize: asNumber(setting(services, "minsize"), "unoise minsize"),
    },
    VSEARCH_TABLE: {
      t: threads,
    },
    marker,
    ITSpart: setting(services, "ITSpart"),
    RDP: {
      memory: `-Xmx${memoryGb}g`,
      custom,
      t: classifierPath,
      c: 0,
      f: "fixrank",
      g: setting(services, "rdp_g"),
    },
    pseudogene_filtering: filterPseudogenes ? "yes" : "no",
    grep_type: grepType,
    taxon1: kept ? `-e ${kept}` : "",
    taxon2: dropped ? `-v ${dropped}` : "",
    removal_type: removalType,
    hmm: "bold.hmm",
    ORFFINDER: {
      g: Number(setting(services, "genetic_code")),
      s: asNumber(setting(services, "orf_start"), "ORF start codon"),
      ml: asNumber(setting(services, "orf_ml"), "minimum ORF length"),
      n: setting(services, "orf_nested") === true ? "true" : "false",
      strand: String(setting(services, "orf_strand")),
    },
    report_type: 2,
  };
}

function assertRdpProperties(database) {
  const fd = fs.openSync(database, "r");
  const buf = Buffer.alloc(8192);
  let head = "";
  try {
    const n = fs.readSync(fd, buf, 0, buf.length, 0);
    head = buf.slice(0, n).toString("utf8");
  } finally {
    fs.closeSync(fd);
  }
  const firstLine = head.split(/\r?\n/).find((line) => line.trim()) || "";
  if (firstLine.startsWith(">") || !head.includes("probabilityList")) {
    throw new MetaWorksConfigError(
      "Select rRNAClassifier.properties from an RDP-trained classifier folder. A FASTA or SINTAX file cannot be used here. Classifier table: https://terrimporter.github.io/MetaWorksSite/#classifier_table",
      "FORM"
    );
  }
}

function prepareMetaWorksConfig({
  inputDir,
  readType,
  services,
  patternOverride,
  memoryBytes,
  cpuCount,
}) {
  if (readType && readType !== "paired_end") {
    throw new MetaWorksConfigError(
      "The MetaWorks ESV workflow is for paired-end reads. Select paired-end when you choose the folder.",
      "FLAT"
    );
  }

  const files = listTopLevelFastqGz(inputDir);
  const readPattern = patternOverride
    ? patternOverride.trim()
    : inferReadPattern(files);
  if (patternOverride) {
    assertPatternMatches(files, readPattern);
  } else if (files.length === 0) {
    inferReadPattern(files);
  }

  const marker = setting(services, "marker");
  let classifierDir = "";
  let classifierContainerPath = "/dev/null";
  if (needsClassifier(marker)) {
    const database = setting(services, "database");
    if (!database || database === "undefined") {
      throw new MetaWorksConfigError(
        "Select the RDP classifier rRNAClassifier.properties file.",
        "FORM"
      );
    }
    if (!fs.existsSync(database)) {
      throw new MetaWorksConfigError(
        `Classifier file not found: ${database}`,
        "FORM"
      );
    }
    assertRdpProperties(database);
    classifierDir = path.dirname(database);
    classifierContainerPath = `/extraFiles/${path.basename(database)}`;
  }

  const config = buildConfigObject(
    services,
    readPattern,
    classifierContainerPath,
    resourcesFromManager(memoryBytes, cpuCount)
  );
  const adapters = buildAdaptersFasta(
    setting(services, "forward_primers"),
    setting(services, "reverse_primers")
  );

  const outputDir = path.join(inputDir, "metaworks_out");
  fs.rmSync(outputDir, { recursive: true, force: true });
  fs.mkdirSync(outputDir, { recursive: true });
  fs.writeFileSync(path.join(outputDir, "adapters.fasta"), adapters);
  fs.writeFileSync(
    path.join(outputDir, "config_ESV.yaml"),
    yaml.dump(config, { lineWidth: -1, noRefs: true, forceQuotes: true })
  );

  return {
    outputDir,
    readPattern,
    classifierDir,
    config,
  };
}

module.exports = {
  MetaWorksConfigError,
  applyPseudogeneDefaults,
  buildAdaptersFasta,
  inferReadPattern,
  listTopLevelFastqGz,
  needsClassifier,
  patternToRegExp,
  prepareMetaWorksConfig,
  reverseComplement,
};
