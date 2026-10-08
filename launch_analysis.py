"""Choose a TIFF and run the existing analysis from a desktop dialog."""
from datetime import datetime
from pathlib import Path
import subprocess
import sys
import tkinter as tk
from tkinter import filedialog, messagebox


def main():
    base = Path(__file__).resolve().parent
    root = tk.Tk()
    root.withdraw()
    filename = filedialog.askopenfilename(
        parent=root, title='분석할 TIFF 파일 선택',
        initialdir=str(base / 'input'),
        filetypes=[('TIFF 이미지', '*.tif *.tiff'), ('모든 파일', '*.*')])
    if not filename:
        root.destroy()
        return
    automatic = messagebox.askyesnocancel(
        '분석 방법 선택',
        '자동으로 블록을 찾을까요?\n\n'
        '예: 자동 검출 (불확실하면 직접 선택)\n'
        '아니요: 직접 네 모서리 선택\n'
        '취소: 종료\n\n'
        '직접 선택 순서: 왼쪽 위 → 오른쪽 위 → 왼쪽 아래 → 오른쪽 아래\n'
        '네 점을 선택한 뒤 Analyze를 누르세요.', parent=root)
    if automatic is None:
        root.destroy()
        return
    output = base / 'results' / (
        Path(filename).stem + '_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    root.destroy()
    result = subprocess.run([
        sys.executable, str(base / 'main.py'), filename,
        '--grid', 'auto' if automatic else 'manual', '--output', str(output)], cwd=base)
    root = tk.Tk()
    root.withdraw()
    if result.returncode == 0:
        messagebox.showinfo('분석 완료', f'결과 저장 위치:\n{output}', parent=root)
    else:
        messagebox.showerror(
            '분석이 완료되지 않았습니다',
            '실행 창의 오류 내용을 확인하세요.\n'
            '선택 창을 닫았다면 분석이 취소된 것입니다.\n\n'
            f'결과가 일부 생성되었다면 저장 위치:\n{output}', parent=root)
    root.destroy()


if __name__ == '__main__':
    main()
