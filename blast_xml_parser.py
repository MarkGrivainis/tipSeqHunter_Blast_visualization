from Bio.Blast.Applications import NcbiblastnCommandline

blastn_cline = NcbiblastnCommandline(task="megablast",
                                     query="temp_sim.fa",
                                     subject="temp_ref.fa", max_target_seqs=100, word_size=28,
                                     penalty=-2,
                                     reward=1,
                                     outfmt=5)
stdout, stderr = blastn_cline()
print(stdout)