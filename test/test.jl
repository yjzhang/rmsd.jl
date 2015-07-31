using RMSD

qui = read_data("../rotation_tests/qui_chr4_fish_nointer.txt")
sen = read_data("../rotation_tests/sen_chr4_fish_nointer.txt")

(_, R, _) = rmsd_rotation(qui, sen)
qui_rot = (R*qui')'
output_txt(qui_rot, "../rotation_tests/qui_rot.txt")
