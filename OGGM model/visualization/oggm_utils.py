import os
import xarray as xr
import geopandas as gpd
from shapely import geometry as shpg
import logging

LOG_FILE = "oggm_aggregate.log" 
# if os.path.exists(LOG_FILE):
#     os.remove(LOG_FILE)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)

def aggregate_to_region(region_info, glacier_gdf, oggm_path_list, config):
    """
    将OGGM模拟结果聚合到指定区域。
    此函数为并行处理的“工作单元”。

    参数:
        region_info (tuple): 一个包含索引和GeoSeries的元组 (index, row)，由 a_gdf.iterrows() 提供。
        glacier_gdf (GeoDataFrame): 冰川边界数据。
        oggm_path_list (list): OGGM数据路径列表。
        config (dict): 包含配置信息的字典，如 'oggm_dir'。
    """
    index, region_row = region_info
    region_boundary = gpd.GeoDataFrame([region_row], crs=glacier_gdf.crs)
    region_name = region_boundary['Name'].values[0]

    try:
        # 筛选区域内的冰川
        # 使用 sjoin 进行空间连接，效率更高
        gdf_sel = gpd.sjoin(glacier_gdf, region_boundary, how="inner", predicate='within') #计算很快
        
        if gdf_sel.empty:
            logging.info(f"【{region_name}】区域内没有冰川，跳过。")
            return f"Skipped: {region_name} (no glaciers)"

        logging.info(f"【{region_name}】区域内冰川数量: {len(gdf_sel)}")
        
        # 为该区域构建输出目录
        output_dir = os.path.join(os.path.dirname(config['oggm_dir']), f'[{index+17}]{region_name}')
        if '/' in output_dir:
            output_dir = output_dir.replace('/', '-')
        os.makedirs(output_dir, exist_ok=True)
        selected_rgi_ids = gdf_sel['RGIId'].tolist()

        for oggm_path in oggm_path_list:
            with xr.open_dataset(oggm_path, chunks='auto') as ds:
                available_rgi_ids = ds['rgi_id'].values # the whole HMA
                valid_rgi_ids = [rid for rid in selected_rgi_ids if rid in available_rgi_ids]
                ds_selected = ds.sel(rgi_id=valid_rgi_ids)
                
                if ds_selected.sizes.get('rgi_id', 0) == 0:
                    logging.warning(f"【{region_name}】在 {os.path.basename(oggm_path)} 中未找到匹配的冰川数据。")
                    continue
                else:
                    logging.info(f"【{region_name}】在 {os.path.basename(oggm_path)} 中找到匹配的冰川数据: {len(valid_rgi_ids)}")
                
                output_file = os.path.join(output_dir, os.path.basename(oggm_path))
                ds_selected.to_netcdf(output_file, compute=True)
                logging.info(f"【{region_name}】成功写入文件: {output_file}")
    
    except Exception as e:
        logging.error(f"处理区域【{region_name}】时发生严重错误: {e}")
        return f"Failed: {region_name} ({e})"

    return f"Success: {region_name}"

