function write_pt(pt)

fileID = fopen('+Sens/pt_params.geo', 'w');
fprintf(fileID, 'tx1_x = %f;\n', pt(1));
fprintf(fileID, 'tx2_x = %f;\n', pt(2));
fprintf(fileID, 'rx1_x = %f;\n', pt(3));
fprintf(fileID, 'rx2_x = %f;\n', pt(4));
fclose(fileID);

end
