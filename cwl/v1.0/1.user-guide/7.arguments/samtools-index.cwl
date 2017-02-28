cwlVersion: v1.0
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: comics/samtools
baseCommand: ["samtools", "index"]
arguments: 
  - valueFrom: '$(runtime.outdir)'
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