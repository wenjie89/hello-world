function [alpha_i, alpha_a] = CalculateAlpha(e, h, option)



    nui_air1=air1(e,h,10);
    nui_morrow=morrowair(e,h,1);

    nua_air1=air1(e,h,2);
    nua_morrow=morrowair(e,h,2);

    mue_air1 = air1(e,h,11);
    mue_morrow = morrowair(e,h,4);


    ve_air1 = -mue_air1.*e;
    ve_morrow = -mue_morrow.*e;


    alpha_i_air1 = nui_air1./abs(ve_air1);
    alpha_a_air1 = nua_air1./abs(ve_air1);

    alpha_i_morrow = nui_morrow./abs(ve_morrow);
    alpha_a_morrow = nua_morrow./abs(ve_morrow);
    
    if strcmp(option,'air1')
        alpha_i = alpha_i_air1;
        alpha_a = alpha_a_air1;
    elseif strcmp(option,'morrowair')
            alpha_i = alpha_i_morrow;
            alpha_a = alpha_a_morrow;
    else
        error('Wrong option!')
    end
       
        

end