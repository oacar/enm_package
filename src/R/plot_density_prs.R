
plot_density_prs <- 			# 
    function
(
 	df , x , xlab= 'Effectiveness', ylab='Density', coord_trans='log10', linesize = 1
)
{
dens_plot <- 
    ggplot()+
    #facet_grid(rows=vars(data))+
#    geom_density()+
    geom_line(data= df%>%filter(group=='real'),
               aes(x={{x}} ,y=..scaled.., color=group),
               stat = 'density')+#, binwidth = 0.1) + 
    geom_line(data= df%>%filter(group=='random'),
               aes(x={{x}} ,y=..scaled..,fill=random_id, color=group),
               stat = 'density',linetype = "dashed" )+#, binwidth = 0.1) + 
#    geom_errorbar(data= df%>%filter(group=='real'), stat = "vline", xintercept = "mean",
#  width=0.8,aes(x = eff_norm , xmax=..x..,xmin=..x..))+
    geom_vline(data= df%>%filter(group=='real'),aes(xintercept = mean({{x}}), color=group),size=linesize, linetype='dotted')+
    geom_vline(data= df%>%filter(group=='random'),aes(xintercept = mean({{x}}), color=group, group=random_id), size=linesize, linetype='dotted')+
    #geom_histogram(binwidth=0.0001)+
    theme_bw()+
    ylab(ylab)+
    xlab(xlab)+
    coord_trans(x = coord_trans)+
scale_x_continuous(
   breaks = scales::trans_breaks("log10", function(x) 10^x),
   labels = scales::trans_format("log10", scales::math_format(10^.x))
 ) +
    scale_color_discrete(labels=c('Rewired\nnetworks','Real\nnetwork'))+
theme(panel.grid.minor = element_blank(), legend.position='bottom', legend.title=element_blank())#+theme_fontsize+
dens_plot
}
#tmp  <- plot_density_prs(pcc_df, eff_norm)
#ggsave('tmp.png', plot=tmp)
