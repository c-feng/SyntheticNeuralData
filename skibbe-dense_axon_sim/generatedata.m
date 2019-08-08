for i = 0:99
    data=mhs_simualte_axions('shape',[300,300,300],'mode','chaos','density',0.5,'max_b',2,'axon_thickness',[4,6],'thick_min',1.5,'flatty_axions',true);
    img = data.Img;
    imgD = data.ImgD;
    gt = data.GT.data;
    con = data.GT.connections;
    % shape_density_maxb_axon_thickness
    filepath = strcat("data/synthetic/chaos_noise/300_0-5_2_4_6/Syn_chaos_300_0-5_2_4_6_", num2str(i, '%03d'), '.mat')
    save(filepath, 'img', 'gt', 'con', 'imgD');
end