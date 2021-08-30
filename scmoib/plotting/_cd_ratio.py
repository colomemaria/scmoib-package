def cd_ratio(df):
    tdf = df.copy()
    tdf['Connected'] = 1 - tdf['cd_ratio']
    tdf.rename(columns={'cd_ratio': 'Disconnected'}, inplace=True)
    tdf[['Disconnected', 'Connected']].plot(kind='bar', stacked=True)