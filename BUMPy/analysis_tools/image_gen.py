import main

template = main.Molecules( 'small_flat_bilayer.pdb')
template.duplicate_laterally(20, 20)
shape = main.shapes.cylinder.gen_shape(template, [0, 0], 50, 200, 0.5, print_intermediates='hal_cyl_pretransform')
shape.write_coordinates('half_cyl_post_transform.pdb')


shape = main.shapes.cylinder.gen_shape(template, [10, 10], 100, 200, 1, print_intermediates='full_cyl_pretransform.pdb')
shape.write_coordinates('full_cyl_post_transform.pdb')

shape = main.shapes.semisphere.gen_shape(template, [10, 10], 100, print_intermediates='semisph_pretransform.pdb')
shape.write_coordinates('semisph_post_transform.pdb')
