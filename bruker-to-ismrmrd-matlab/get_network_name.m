function [ str_network_imagerie, str_network_perso ] = get_network_name( str_user )
%fonction qui permet de travailler sur les odirnateurs ubuntu et le
%calculateur

if (strcmp(str_user,'valery'))
    
    str_network_imagerie='/Reseau/Imagerie/';
    str_network_perso='/Reseau/Valery/';
    
elseif (strcmp(str_user,'valeryo'))
    
   str_network_imagerie='/mount/Imagerie/';
    str_network_perso='/mount/valery.ozenne/';    
    
elseif (strcmp(str_user,'juliem'))
    
   str_network_imagerie='/mount/Imagerie/';
    str_network_perso='/mount/julie.magat/';
   
elseif (strcmp(str_user,'kylianh'))

str_network_imagerie='/mount/Imagerie/';
str_network_perso='/mount/kylian.haliot/';
    
else
   
    str_network_imagerie='/Reseau/Imagerie/';
    str_network_perso='/Reseau/Valery/';
    
end

end

