"use strict";
(() => {
  const search = document.querySelector("#project-search");
  const cards = [...document.querySelectorAll(".project-card")];
  const filters = [...document.querySelectorAll("[data-filter]")];
  let domain = "all";
  function applyFilter() {
    const words = (search?.value || "").trim().toLocaleLowerCase().split(/\s+/).filter(Boolean);
    let visible = 0;
    for (const card of cards) {
      const content = card.dataset.search.toLocaleLowerCase();
      const matches = (domain === "all" || card.dataset.domain === domain) && words.every(word => content.includes(word));
      card.hidden = !matches;
      visible += Number(matches);
    }
    const count = document.querySelector("#result-count");
    const empty = document.querySelector("#empty-state");
    if (count) count.textContent = `共 ${visible} 项成果`;
    if (empty) empty.hidden = visible !== 0;
    for (const button of filters) button.setAttribute("aria-pressed", String(button.dataset.filter === domain));
  }
  search?.addEventListener("input", applyFilter);
  for (const button of filters) button.addEventListener("click", () => { domain = button.dataset.filter; applyFilter(); });
  document.querySelector("#clear-search")?.addEventListener("click", () => { search.value = ""; domain = "all"; applyFilter(); search.focus(); });
  for (const button of document.querySelectorAll("[data-copy]")) {
    button.addEventListener("click", async () => {
      const code = button.parentElement.querySelector("pre code");
      const status = button.parentElement.querySelector(".copy-status");
      try {
        await navigator.clipboard.writeText(code.textContent);
        status.textContent = "命令已复制；请在对应项目目录执行。";
      } catch {
        const range = document.createRange();
        range.selectNodeContents(code);
        const selection = window.getSelection();
        selection.removeAllRanges();
        selection.addRange(range);
        status.textContent = "已选中命令，可按 Ctrl+C 复制。";
      }
    });
  }
})();
