cwlVersion: v1.0
class: CommandLineTool
requirements:
  - class: InlineJavascriptRequirement
hints:
  DockerRequirement:
    dockerPull: comics/samtools
baseCommand: ["samtools", "index"]
arguments: 
  - valueFrom: $(runtime.outdir + "/" + inputs.bam.path.split('/').slice(-1) + ".bai")
    position: 2
inputs:
  bam:
    type: File
    inputBinding:
      position: 1
outputs:
  bai:
    type: File
    outputBinding:
      glob: "*.bai"