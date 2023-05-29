`GATK best practices workflow`_ Pipeline summary

=============================================
Reference
=============================================
1. Reference genome related files and GTAK budnle files (GATK_)
2. VEP Variarition annotation files (VEP_)

=============================================
Prepare
=============================================
1. Adapter trimming (Fastp_)
2. Aligner (`BWA mem2`_)
3. Mark duplicates (samblaster_)
4. Generates recalibration table for Base Quality Score Recalibration (BaseRecalibrator_)
5. Apply base quality score recalibration (ApplyBQSR_)

=============================================
Quality control report
=============================================
1. Fastp report (MultiQC_)
2. Alignment report (MultiQC_)

=============================================
Somatic short variants (SNV+INDEL)
=============================================
1. Call somatic SNVs and indels via local assembly of haplotypes (Mutect2_)
2. Tabulates pileup metrics for inferring contamination (GetPileupSummaries_)
3. Calculate the fraction of reads coming from cross-sample contamination (CalculateContamination_)
4. Get the maximum likelihood estimates of artifact prior probabilities in the orientation bias mixture model filter (LearnReadOrientationModel_)
5. Filter somatic SNVs and indels called by Mutect2 (FilterMutectCalls_)
6. Merge all the VCF files (Picard_)

=============================================
Annotation
=============================================
Annotate variant calls with VEP (VEP_)

.. _GATK best practices workflow: https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows
.. _GATK: https://software.broadinstitute.org/gatk/
.. _VEP: https://www.ensembl.org/info/docs/tools/vep/index.html
.. _fastp: https://github.com/OpenGene/fastp
.. _BWA mem2: http://bio-bwa.sourceforge.net/
.. _samblaster: https://github.com/GregoryFaust/samblaster
.. _BaseRecalibrator: https://gatk.broadinstitute.org/hc/en-us/articles/13832708374939-BaseRecalibrator
.. _ApplyBQSR: https://github.com/GregoryFaust/samblaster
.. _Picard: https://broadinstitute.github.io/picard
.. _Mutect2: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-Mutect2
.. _GetPileupSummaries: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-GetPileupSummaries
.. _CalculateContamination: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-CalculateContamination
.. _LearnReadOrientationModel: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-LearnReadOrientationModel
.. _FilterMutectCalls: https://gatk.broadinstitute.org/hc/en-us/articles/13832694334235-FilterMutectCalls
.. _MultiQC: https://multiqc.info
