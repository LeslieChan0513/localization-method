for p = 1:Itermax_a
    val_sub = 0;
    %% Recording
    for i = 1:anchors_n
        val_sub = val_sub + 1/2 * ((c * (ti(i,1) - t0_kk) - norm(all_nodes.all(:,i) - y_k))^2);
    end
    %% ��ѭ����������
    grad = zeros(2,1);
    for i = 1:anchors_n
        grad = grad + (c * (ti(i,1) - t0_kk) - norm(all_nodes.all(:,i) - y_k)) * (all_nodes.all(:,i) - y_k)/(norm(all_nodes.all(:,i) - y_k));
    end
    if (norm(grad) < eps_p)
         break;
    end
    %% ��������
    d = -grad;
    %% ������ʲ���
    yn = 0;
    alpha_a = 1;                                                        %��ʼ����
    yn = y_k + alpha_a * d;
    fn_sub = 0;
    for i = 1:anchors_n
        fn_sub =  fn_sub + 1/2 *((c * (ti(i,1) - t0_kk) - norm(all_nodes.all(:,i) - yn))^2);
    end
    while (fn_sub > val_sub + alpha_a * sigma * grad' * d )                   %Armijo��������
        alpha_a = alpha_a * 0.618;                                    %��0.618����
        yn = y_k + alpha_a * d;                                       %���²�������ֵ�ж��Ƿ������������
        fn_sub = 0;
        for i = 1:anchors_n
            fn_sub =  fn_sub + 1/2 *((c * (ti(i,1) - t0_kk) - norm(all_nodes.all(:,i) - yn))^2);
        end
    end
    %% ����yֵ         
     y_kk =  y_k + alpha_a * d;                                           %ʹ�õõ����²�������ֵ
     break
end
y_bb = y_kk;
y_ini = y_k;
y_new = y_kk;