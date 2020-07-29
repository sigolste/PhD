function txt = progress_bar(var,var_tot,var_name)

L = 50;

per_num = floor(var/var_tot * L);

progress = repmat('=',1,per_num);
spaces = repmat(' ',1,L-per_num);

txt = [var_name,'\n[',progress,spaces,'] ',num2str(var),'/',num2str(var_tot),'\n'];

end